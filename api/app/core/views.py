import csv
from datetime import datetime
from uuid import uuid4
import re
import math
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
    WellDetail,
    PlateDetail,
    MeasurementFeature,
)
from .serializers import (
    PlateSerializer,
    WellSerializer,
    PlateMappingSerializer,
    SimpleExperimentSerializer,
    ExperimentSerializer,
    ProjectSerializer,
    SimplePlateTemplateSerializer,
    WellDetailSerializer,
    PlateDetailSerializer,
)

from django.views.generic import View
from django.conf import settings


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
        expression = request.data.get("expression")
        new_label = request.data.get("new_label")
        now = datetime.now().replace(microsecond=0)

        current_plate = self.get_object()

        # find indices of corresponding measurements in the first well of the current plate
        # we need that because the timestamps are not the same across the experiment plate
        # so that we will apply the formula to other plates according to the order of the
        # measurements passed in the formula from the ui
        current_plate_first_well = current_plate.wells.first()
        first_well_measurements_indices = {}
        for index, measurement in enumerate(
            current_plate_first_well.measurements.all()
        ):
            measurement_combined_label = (
                f""
                f"{measurement.label}"
                f"_{measurement.measurement_timestamp.isoformat().split('+')[0]}"
            )
            first_well_measurements_indices[measurement_combined_label] = index

        experiment = current_plate.experiment
        experiment_plates = experiment.plates.all()
        measurement_feature, _ = MeasurementFeature.objects.get_or_create(
            abbrev=new_label
        )

        try:
            for plate in experiment_plates:
                wells = plate.wells.all()
                for well in wells:
                    well_measurements = []
                    for measurement in well.measurements.all():
                        well_measurements.append(measurement.value)
                    new_formula = expression.replace("ln(", "math.log(")

                    for key, value in first_well_measurements_indices.items():
                        new_formula = new_formula.replace(
                            key, str(well_measurements[value])
                        )

                    try:
                        result = eval(new_formula)
                        if not result:
                            # if the result is None or 0, we will set it to 0, for example if the
                            # value = 0, it's log can be -inf, so we will set it to 0
                            log.error(
                                f"Result is None or 0. Setting result to 0. Formula: {new_formula}"
                            )
                            result = 0
                    except ZeroDivisionError:
                        log.error("Division by zero occurred. Setting result to 0")
                        result = 0

                    Measurement.objects.create(
                        well=well,
                        label=new_label,
                        value=result,
                        measurement_timestamp=now,
                        identifier="",
                        feature=measurement_feature,
                    )
        except (ValueError, TypeError, NameError) as e:
            return Response({"error": str(e)}, status.HTTP_400_BAD_REQUEST)
        else:
            return Response(status.HTTP_200_OK)

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)


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


class WellDetailViewSet(viewsets.ModelViewSet):
    serializer_class = WellDetailSerializer
    queryset = WellDetail.objects.all()


class PlateDetailViewSet(viewsets.ModelViewSet):
    serializer_class = PlateDetailSerializer
    queryset = PlateDetail.objects.all()
