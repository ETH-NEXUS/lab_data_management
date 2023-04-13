import csv
from datetime import datetime

from uuid import uuid4
import re
from friendlylog import colored_logger as log

from os import environ
import os

from compoundlib.serializers import SimpleCompoundLibrarySerializer
from django.core.exceptions import ValidationError
from django.db import IntegrityError
from django.db.models import Prefetch, Q
from django.http import Http404
from rest_framework import viewsets, views, status
from rest_framework.decorators import action
from rest_framework.response import Response
from django.http import HttpResponse
import mimetypes
from django.utils.translation import gettext as _

from .models import (
    Well,
    Plate,
    Measurement,
    WellWithdrawal,
    WellCompound,
    PlateMapping,
    Experiment,
    BarcodeSpecification,
    PlateDimension,
    Project,
    MeasurementFeature,
    PlateDetail,
    WellDetail,
    ExperimentDetail,
)
from .serializers import (
    PlateSerializer,
    WellSerializer,
    PlateMappingSerializer,
    SimpleExperimentSerializer,
    ExperimentSerializer,
    ProjectSerializer,
    SimplePlateTemplateSerializer,
)

from django.views.generic import View
from django.conf import settings


def mean_time_point(dt_strings):
    dt_array = [datetime.fromisoformat(dt_str) for dt_str in dt_strings]
    timestamps = [dt.timestamp() for dt in dt_array]
    avg_timestamp = sum(timestamps) / len(timestamps)
    avg_datetime = datetime.fromtimestamp(avg_timestamp)
    return avg_datetime


class ProjectViewSet(viewsets.ModelViewSet):
    serializer_class = ProjectSerializer

    def get_queryset(self):
        experiments = Prefetch("experiments", queryset=Experiment.objects.all())
        plates = Prefetch("plates", queryset=Plate.objects.all().order_by("barcode"))
        return (
            Project.objects.all().prefetch_related(experiments).prefetch_related(plates)
        )


class PlateViewSet(viewsets.ModelViewSet):
    def get_serializer_class(self):
        # if self.action == 'list':
        #     return PlateListSerializer
        # else:
        return PlateSerializer

    def get_queryset(self):
        measurements = Prefetch("measurements", queryset=Measurement.objects.all())
        withdrawals = Prefetch(
            "withdrawals",
            queryset=WellWithdrawal.objects.select_related("target_well").all(),
        )
        donors = Prefetch(
            "donors",
            queryset=WellWithdrawal.objects.select_related("well").all(),
        )
        well_compounds = Prefetch(
            "well_compounds",
            queryset=WellCompound.objects.select_related(
                "compound", "compound__library"
            ).all(),
        )
        wells = Prefetch(
            "wells",
            queryset=Well.objects.select_related("sample", "type")
            .order_by("position")
            .prefetch_related(well_compounds)
            .prefetch_related(withdrawals)
            .prefetch_related(donors)
            .prefetch_related(measurements),
        )
        return Plate.objects.select_related(
            "dimension", "experiment", "library", "template"
        ).prefetch_related(wells)

    @action(detail=False, methods=["get"])
    def barcodes(self, request):
        """Returns an array of barcodes"""
        library = request.GET.get("library")
        experiment = request.GET.get("experiment")
        template = request.GET.get("template")
        predicate = Q()
        if library:
            predicate |= Q(library__isnull=(library.lower() != "true"))
        if experiment:
            predicate |= Q(experiment__isnull=(experiment.lower() != "true"))
        if template:
            predicate |= Q(template__isnull=(template.lower() != "true"))
        return Response(
            [
                {
                    "label": plate.barcode
                    if plate.template is None
                    else " / ".join(plate.barcode.replace("__TEMPL__", "").split("_")),
                    "value": plate.id,
                    "library": SimpleCompoundLibrarySerializer(plate.library).data
                    if plate.library
                    else None,
                    "experiment": SimpleExperimentSerializer(plate.experiment).data
                    if plate.experiment
                    else None,
                    "template": SimplePlateTemplateSerializer(plate.template).data
                    if plate.template
                    else None,
                }
                for plate in Plate.objects.filter(predicate)
            ]
        )

    @action(detail=True, methods=["post"])
    def apply_template(self, request, pk=None):
        """Applies a template plate"""
        plate = self.get_object()
        template_plate_id = request.data.get("template")
        if template_plate_id:
            template_plate = Plate.objects.get(pk=template_plate_id)
            plate = plate.apply_template(template_plate)
            return Response(PlateSerializer(plate).data, status.HTTP_200_OK)
        else:
            raise Http404("Parameter 'template' is required.")

    @action(detail=True, methods=["post"])
    def add_new_measurement(self, request, pk=None):
        import math

        used_labels = request.data.get("used_labels")
        new_label = request.data.get("new_label")
        expression = request.data.get("expression").replace("ln(", "math.log(")

        plate_id = request.data.get("plate_id")
        experiment_id = request.data.get("experiment_id")
        separate_time_series_points = request.data.get("separate_time_series_points")

        if plate_id:
            current_plate = Plate.objects.get(id=plate_id)
            current_plate_details = PlateDetail.objects.get(id=plate_id)

            self.__add_new_measurement_to_plate(
                current_plate,
                current_plate_details,
                used_labels,
                new_label,
                expression,
                separate_time_series_points,
            )
        elif experiment_id:
            experiment = Experiment.objects.get(id=experiment_id)
            plates = Plate.objects.filter(experiment=experiment)
            for plate in plates:
                current_plate_details = PlateDetail.objects.get(id=plate.id)
                self.__add_new_measurement_to_plate(
                    plate, current_plate_details, used_labels, new_label, expression
                )

        return Response(status.HTTP_200_OK)

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)

    def __evaluate_expression(self, new_expression):
        import math

        try:
            result = eval(new_expression)
            if not result:
                log.warning(
                    f"Result is None or 0. Setting result to 0. Formula: "
                    f"{new_expression}"
                )
                result = 0
        except ZeroDivisionError:
            log.error("Division by zero occurred. Setting result to 0")
            result = 0
        return result

    def __add_new_measurement_to_plate(
        self,
        current_plate,
        current_plate_details,
        used_labels,
        new_label,
        expression,
        separate_time_series_points,
    ):
        new_measurement_timestamp, time_series_support = self.__create_new_timestamp(
            current_plate_details
        )
        wells = current_plate.wells.all()
        measurement_feature, _ = MeasurementFeature.objects.get_or_create(
            abbrev=new_label
        )
        if separate_time_series_points:
            self.__separate_time_series_points(
                used_labels,
                new_label,
                expression,
                wells,
                measurement_feature,
                new_measurement_timestamp,
            )
        else:
            self.__all_time_series_points(
                used_labels,
                new_label,
                expression,
                wells,
                measurement_feature,
                time_series_support,
                new_measurement_timestamp,
            )

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        log.info(
            f"New measurement {new_label} added to plate {current_plate.id} with barcode {current_plate.barcode}"
        )

    def __separate_time_series_points(
        self,
        used_labels,
        new_label,
        expression,
        wells,
        measurement_feature,
        new_measurement_timestamp,
    ):
        measurement_objects = []
        for item in used_labels:
            label, timestamp = item.split("-->")
            label = label.strip().lstrip()
            timestamp = timestamp.strip().lstrip()
            measurement_objects.append(
                {"label": label, "timestamp": timestamp, "combined_label": item}
            )

        for well in wells:
            well_measurements = well.measurements.all()
            new_expression = expression
            for measurement in well_measurements:
                for measurement_object in measurement_objects:
                    if measurement.label == measurement_object[
                        "label"
                    ] and measurement_object["timestamp"].split("+")[
                        0
                    ] == measurement.measured_at.strftime(
                        "%Y-%m-%dT%H:%M:%S%z"
                    ):
                        new_expression = new_expression.replace(
                            measurement_object["combined_label"],
                            str(measurement.value),
                        )

            measurement, _ = self.__create_measurement(
                well,
                new_label,
                new_expression,
                new_measurement_timestamp,
                "",
                measurement_feature,
            )

    def __all_time_series_points(
        self,
        used_labels,
        new_label,
        expression,
        wells,
        measurement_feature,
        time_series_support,
        new_measurement_timestamp,
    ):
        for well in wells:
            well_measurements = well.measurements.all()
            if time_series_support:
                for measurement in well_measurements:
                    current_measurement_time = measurement.measured_at
                    same_time_measurements_data = {}
                    for _measurement in well_measurements:
                        if _measurement.measured_at == current_measurement_time:
                            same_time_measurements_data[
                                _measurement.label
                            ] = _measurement.value

                    new_expression = expression
                    if (
                        len(same_time_measurements_data) > 0
                        and measurement.label in used_labels
                    ):
                        for key in same_time_measurements_data.keys():
                            new_expression = new_expression.replace(
                                key, str(same_time_measurements_data[key])
                            )

                        self.__create_measurement(
                            well,
                            new_label,
                            new_expression,
                            current_measurement_time,
                            "",
                            measurement_feature,
                        )

            else:
                new_expression = expression
                for measurement in well_measurements:
                    new_expression = new_expression.replace(
                        measurement.label, str(measurement.value)
                    )

                self.__create_measurement(
                    well,
                    new_label,
                    new_expression,
                    new_measurement_timestamp,
                    "new_value",
                    measurement_feature,
                )

    def __create_new_timestamp(self, current_plate_details):
        now = datetime.now().replace(microsecond=0)
        time_series_support = False
        new_measurement_timestamp = now

        if len(current_plate_details.measurement_labels) == 1:
            measurement_label = current_plate_details.measurement_labels[0]
            if len(current_plate_details.measurement_timestamps[measurement_label]) > 1:
                time_series_support = True
            else:
                new_measurement_timestamp = (
                    current_plate_details.measurement_timestamps[measurement_label][0]
                )
        else:
            measurement_label = list(
                current_plate_details.measurement_timestamps.keys()
            )[0]
            if len(current_plate_details.measurement_timestamps[measurement_label]) > 1:
                time_series_support = True
            else:
                timestamps = []
                for key in current_plate_details.measurement_timestamps.keys():
                    timestamps.append(
                        current_plate_details.measurement_timestamps[key][0]
                    )
                new_measurement_timestamp = mean_time_point(timestamps)
        return new_measurement_timestamp, time_series_support

    def __create_measurement(
        self, well, label, new_expression, measured_at, identifier, feature
    ):
        value = self.__evaluate_expression(new_expression)
        measurement, _ = Measurement.objects.get_or_create(
            well=well,
            label=label,
            measured_at=measured_at,
            defaults={"value": value, "identifier": identifier, "feature": feature},
        )
        return measurement, _


class WellViewSet(viewsets.ModelViewSet):
    serializer_class = WellSerializer
    queryset = Well.objects.all()

    @action(detail=True, methods=["get"])
    def chain(self, request, pk=None):
        def node_key(well):
            return f"{well.plate.barcode}_{well.hr_position}"

        def node_name(well):
            return f"{well.plate.barcode}: {well.hr_position}"

        def appendDonors(nodes, edges, well):
            for donor in well.donors.all():
                nodes.update({node_key(donor.well): {"name": node_name(donor.well)}})
                edge_key = str(uuid4())
                edges.update(
                    {
                        edge_key: {
                            "source": node_key(donor.well),
                            "target": node_key(well),
                            "label": donor.amount,
                        }
                    }
                )
                appendDonors(nodes, edges, donor.well)

        def appendWithdrawals(nodes, edges, well):
            for withdrawal in well.withdrawals.all():
                if withdrawal.target_well:
                    nodes.update(
                        {
                            node_key(withdrawal.target_well): {
                                "name": node_name(withdrawal.target_well)
                            }
                        }
                    )
                    edge_key = str(uuid4())
                    edges.update(
                        {
                            edge_key: {
                                "source": node_key(well),
                                "target": node_key(withdrawal.target_well),
                                "label": withdrawal.amount,
                            }
                        }
                    )
                    appendWithdrawals(nodes, edges, withdrawal.target_well)

        well = Well.objects.get(pk=pk)
        nodes = {node_key(well): {"name": node_name(well), "root": True}}
        edges = {}
        appendDonors(nodes, edges, well)
        appendWithdrawals(nodes, edges, well)

        return Response({"nodes": nodes, "edges": edges}, status=status.HTTP_200_OK)


class PlateMappingViewSet(viewsets.ModelViewSet):
    serializer_class = PlateMappingSerializer
    queryset = PlateMapping.objects.all()


class MappingPreviewView(views.APIView):
    def post(self, request, format=None):
        delimiter = request.GET.get("delimiter") or ","
        quotechar = request.GET.get("quotechar") or '"'
        for _file in request.data:
            with open(_file, "r") as file:
                reader = csv.DictReader(file, delimiter=delimiter, quotechar=quotechar)
                data = [line for line in reader]
            # We expect only one file!
            break

        return Response(data, status=status.HTTP_200_OK)


class ExperimentViewSet(viewsets.ModelViewSet):
    serializer_class = ExperimentSerializer
    queryset = Experiment.objects.all()

    @action(detail=False, methods=["post"])
    def barcodes(self, request):
        """Saves a barcode specifications for an experiment"""
        data = request.data
        experiment = Experiment.objects.get(pk=data["experiment_id"])
        barcode_specification = BarcodeSpecification(
            prefix=data["prefix"],
            number_of_plates=data["number_of_plates"],
            sides=data["sides"],
            experiment=experiment,
        )
        barcode_specification.save()
        return Response(status=status.HTTP_200_OK)

    @action(detail=False, methods=["post"])
    def bulk_apply_template(self, request, pk=None):
        """Applies a template plate to all the plate of the experiment"""
        experiment = Experiment.objects.get(pk=request.data.get("experiment_id"))
        template_plate_id = request.data.get("template")
        if template_plate_id:
            template_plate = Plate.objects.get(pk=template_plate_id)
            plates = Plate.objects.filter(experiment=experiment)
            for plate in plates:
                plate.apply_template(template_plate)
            return Response(status.HTTP_200_OK)
        else:
            raise Http404(
                _("Parameters 'template' and 'experiment_id' are " "required.")
            )

    @action(detail=False, methods=["post"])
    def bulk_add_plates(self, request):
        try:
            experiment_id = int(request.data["experiment_id"])
            barcode_specification_id = int(request.data["barcode_specification_id"])
            plate_dimension_id = int(request.data["plate_dimension_id"])
            experiment = Experiment.objects.get(pk=experiment_id)
            barcode_specification = BarcodeSpecification.objects.get(
                pk=barcode_specification_id
            )
            plate_dimension = PlateDimension.objects.get(pk=plate_dimension_id)
            number_of_plates = barcode_specification.number_of_plates
        except Experiment.DoesNotExist:
            return Response(
                {"error": _("Experiment not found")},
                status=status.HTTP_404_NOT_FOUND,
            )
        except BarcodeSpecification.DoesNotExist:
            return Response(
                {"error": _("Barcode specification not found")},
                status=status.HTTP_404_NOT_FOUND,
            )
        except PlateDimension.DoesNotExist:
            return Response(
                {"error": _("Plate dimension not found")},
                status=status.HTTP_404_NOT_FOUND,
            )

        plates = []
        for i in range(number_of_plates):
            plate = Plate(
                barcode=barcode_specification.get_barcode_by_number(i + 1),
                experiment=experiment,
                dimension=plate_dimension,
            )
            plates.append(plate)

        try:
            Plate.objects.bulk_create(plates)
        except IntegrityError:
            return Response(
                {
                    "error": _(
                        "Could not save plates to database. Probably you have already added the plates "
                        "with these barcode prefix to the current experiment"
                    )
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
        except ValidationError:
            return Response(
                {"error": _("Invalid plate data")},
                status=status.HTTP_400_BAD_REQUEST,
            )

        return Response(
            {"success": _("Plates added successfully")}, status=status.HTTP_200_OK
        )


class VersionView(views.APIView):
    def get(self, request, format=None):
        return Response({"version": environ.get("GIT_VERSION", "N/A")})


class DocsView(View):
    def get(self, request, uri, **kwargs):
        docs_dir = os.path.join(settings.BASE_DIR, "docs", "site")
        if uri == "":
            uri = "index.html"
        file_path = os.path.join(docs_dir, uri)

        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "index.html")
        if not os.path.isfile(file_path):
            raise Http404("File not found")

        with open(file_path, "r", encoding="utf8", errors="replace") as f:
            content = f.read()

        # Replace all emojis with an empty string
        content = re.sub(r"[^\x00-\x7F]+", "", content)
        mime_type = mimetypes.guess_type(file_path)
        return HttpResponse(content, content_type=mime_type[0])
