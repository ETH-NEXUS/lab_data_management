import csv
import json
import mimetypes
import os
import re
from datetime import datetime
from os import environ
from uuid import uuid4
import subprocess
import traceback

from django.conf import settings
from django.contrib.auth import authenticate, login, logout
from django.core.exceptions import ValidationError
from django.core.handlers.wsgi import WSGIRequest
from django.db import IntegrityError
from django.db.models import Prefetch, Q
from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from django.utils.translation import gettext as _
from django.views.decorators.csrf import ensure_csrf_cookie, csrf_exempt
from django.views.generic import View
from rest_framework import viewsets, views, status
from rest_framework.decorators import action
from django.http import JsonResponse
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from ldm.ldm import get_experiment_measurements


from compoundlib.serializers import SimpleCompoundLibrarySerializer
from helpers.logger import logger
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
    PlateInfo,
)
from .serializers import (
    PlateSerializer,
    WellSerializer,
    PlateMappingSerializer,
    SimpleExperimentSerializer,
    ExperimentSerializer,
    ProjectSerializer,
    SimplePlateTemplateSerializer,
    ExperimentDetail,
)

GLOBAL_NOW = datetime.now().replace(microsecond=0)


class CsrfCookieView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request: WSGIRequest, *args, **kwargs):
        return JsonResponse({"details": _("CSRF cookie set")})


class LoginView(View):
    def post(self, request: WSGIRequest, *args, **kwargs):
        data = json.loads(request.body)
        username = data.get("username")
        password = data.get("password")

        if username is None or password is None:
            return JsonResponse(
                {"detail": _("Please provide username and password.")}, status=400
            )

        user = authenticate(username=username, password=password)

        if user is None:
            return JsonResponse({"detail": _("Invalid credentials.")}, status=400)

        login(request, user)
        return JsonResponse({"detail": _("Successfully logged in.")})


class LogoutView(View):
    def get(self, request: WSGIRequest, *args, **kwargs):
        if not request.user.is_authenticated:
            return JsonResponse({"detail": _("You're not logged in.")}, status=400)

        logout(request)
        return JsonResponse({"detail": _("Successfully logged out.")})


def mean_time_point(dt_strings):
    try:
        dt_array = [datetime.fromisoformat(dt_str) for dt_str in dt_strings]
        timestamps = [dt.timestamp() for dt in dt_array]
        avg_timestamp = sum(timestamps) / len(timestamps)
        avg_datetime = datetime.fromtimestamp(avg_timestamp)
        return avg_datetime
    except Exception as e:
        logger.error(f"Error calculating mean time point: {e}")
        traceback.print_exc()
        return GLOBAL_NOW


class ProjectViewSet(viewsets.ModelViewSet):
    serializer_class = ProjectSerializer

    def get_queryset(self):
        experiments = Prefetch("experiments", queryset=Experiment.objects.all())
        plates = Prefetch("plates", queryset=Plate.objects.all().order_by("barcode"))
        return (
            Project.objects.all().prefetch_related(experiments).prefetch_related(plates)
        )


@api_view(["POST"])
def add_control_layout(request):
    """Add a selected control layout to a project"""

    data = request.data
    project_id = data.get("project_id")
    barcode_new = data.get("barcode_new")
    barcode_old = data.get("barcode_old")

    try:
        project = Project.objects.get(id=project_id)
        plate_old = Plate.objects.get(barcode=barcode_old)

        plate_new = Plate(
            barcode=barcode_new,
            dimension=plate_old.dimension,
            library=plate_old.library,
            template=plate_old.template,
            is_control_plate=True,
            project=project,
        )
        plate_new.save()
        plate_old.copy(plate_new, map_type=True)
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        return Response(status=status.HTTP_200_OK)

    except Project.DoesNotExist:
        logger.error(f"Project not found: {project_id}")
        return Response(
            {"error": "Project not found"}, status=status.HTTP_404_NOT_FOUND
        )
    except Plate.DoesNotExist:
        logger.error(f"Old plate not found: {barcode_old}")
        return Response(
            {"error": "Old plate not found"}, status=status.HTTP_404_NOT_FOUND
        )
    except Exception as e:
        traceback.print_exc()
        logger.error(f"Error adding control layout: {e}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


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
            queryset=WellCompound.objects.select_related("compound").all(),
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
        if request.GET.get("barcode"):
            plate = Plate.objects.get(barcode=request.GET.get("barcode"))
            experiment = plate.experiment
            if experiment:
                project = experiment.project
                project_plates = Plate.objects.filter(project=project)

                return Response(
                    [
                        {"label": plate.barcode, "value": plate.id}
                        for plate in project_plates
                    ]
                )
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

    from rest_framework.exceptions import ValidationError

    @action(detail=True, methods=["post"])
    def apply_template(self, request, pk=None):
        """Applies a template plate"""
        apply_to_all_experiment_plates = request.data.get(
            "apply_to_all_experiment_plates"
        )
        template_plate_id = request.data.get("template")
        print("_______________________________________________")
        print(f"apply_to_all_experiment_plates: {apply_to_all_experiment_plates}")
        print(f"template_plate_id: {template_plate_id}")
        if template_plate_id is None:
            raise ValidationError(
                {"template": ["This field is required."]}, code="invalid"
            )

        template_plate = Plate.objects.get(pk=template_plate_id)
        plate = self.get_object()

        if apply_to_all_experiment_plates:
            plates = Plate.objects.filter(experiment=plate.experiment)
            for _plate in plates:
                _plate.apply_template(template_plate)
        else:
            plate.apply_template(template_plate)

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

        return Response(PlateSerializer(plate).data, status=status.HTTP_200_OK)

    @action(detail=True, methods=["post"])
    def add_new_measurement(self, request, pk=None):

        used_labels = request.data.get("used_labels")
        new_label = request.data.get("new_label")
        expression = request.data.get("expression").replace("log(", "math.log10(")

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
                    plate,
                    current_plate_details,
                    used_labels,
                    new_label,
                    expression,
                    separate_time_series_points,
                )
        ExperimentDetail.refresh(concurrently=True)
        return Response(status.HTTP_200_OK)

    def filter_queryset(self, queryset):
        return super().filter_queryset(queryset)

    def __evaluate_expression(self, new_expression):

        try:
            result = eval(new_expression)
            if not result:
                logger.warning(
                    f"Result is None or 0. Setting result to 0. Formula: "
                    f"{new_expression}"
                )
                result = 0
        except ZeroDivisionError:
            logger.critical("Division by zero occurred. Setting result to 0")
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
        logger.info(
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

    def destroy(self, request, *args, **kwargs):
        well = self.get_object()
        well.delete()

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

        return Response(status=status.HTTP_200_OK)

    @action(detail=True, methods=["get"])
    def mark_as_invalid(self, request, pk=None):
        well = self.get_object()
        well.is_invalid = True
        well.save()
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)
        return Response(status=status.HTTP_200_OK)

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
                logger.debug(f"Applying template to plate {plate.barcode}")
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


@csrf_exempt
def generate_pdf_report(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            notebook_path = data.get("notebook_path")
            experiment = data.get("experiment")
            label = data.get("label")
            if not notebook_path:

                notebook_path = "/notebooks/input/general.ipynb"

            cmd = [
                "python",
                "/root/.ipython/profile_default/report_generator.py",
                "--notebook_path",
                notebook_path,
                "--experiment",
                experiment,
                "--label",
                label,
            ]
            subprocess.run(cmd, check=True)
            return JsonResponse({"status": "Report generated successfully"})
        else:
            return JsonResponse({"error": "Invalid request method"}, status=400)
    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
@api_view(["POST"])
def list_files(request):
    try:
        experiment = request.data.get("experiment")
        notebooks_dir = request.data.get("notebooks_dir")
        file_format = request.data.get("file_format")
        print(experiment, notebooks_dir, file_format)
        if not notebooks_dir:
            notebooks_dir = f"/notebooks/output/{experiment}"
        if not file_format:
            file_format = ".pdf"
        notebooks = []
        if os.path.exists(notebooks_dir):
            for file in os.listdir(notebooks_dir):
                if file.endswith(file_format):
                    full_path = os.path.join(notebooks_dir, file)
                    notebooks.append(full_path)
        return Response({"notebooks": notebooks}, status=status.HTTP_200_OK)
    except Exception as e:
        print(e)
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


@csrf_exempt
def download_pdf_report(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            path = data.get("path")

            if not path:
                return JsonResponse({"error": "Path not provided"}, status=400)
            if not os.path.exists(path):
                return JsonResponse({"error": "File not found"}, status=404)
            with open(path, "rb") as f:
                response = HttpResponse(f, content_type="application/pdf")
                response[
                    "Content-Disposition"
                ] = f'attachment; filename="{os.path.basename(path)}"'
                return response
        else:
            return JsonResponse({"error": "Invalid request method"}, status=400)
    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
def download_csv_data(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            label = data.get("label")
            experiment = data.get("experiment")
            df = get_experiment_measurements(experiment, label)
            response = HttpResponse(content_type="text/csv")
            response["Content-Disposition"] = f'attachment; filename="{label}.csv"'
            df.to_csv(response, index=False)
            return response

    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)


def get_existing_plate_infos(experiment_id):
    plate_info = []
    plate_infos = PlateInfo.objects.filter(experiment=experiment_id)
    if plate_infos:
        for item in plate_infos:
            obj = {
                "plate_barcode": item.plate.barcode,
                "lib_plate_barcode": item.lib_plate_barcode,
                "measurement_label": item.label,
                "replicate": item.replicate,
                "measurement_timestamp": item.measurement_time,
                "cell_type": item.cell_type,
                "condition": item.condition,
            }
            plate_info.append(obj)
    return plate_info


def find_withdrawal_well(start_index, wells, total_columns):
    """
    If the middle well is empty, we look for the closest well with withdrawals up and down.
    """
    if len(WellWithdrawal.objects.filter(target_well=wells[start_index])) > 0:
        logger.debug(f"Number of wells: {len(wells)}")
        logger.debug(f"Starting index: {start_index}")
        return wells[start_index], WellWithdrawal.objects.filter(
            target_well=wells[start_index]
        )
    for offset in range(total_columns):
        indices_to_check = [start_index + i * total_columns for i in range(-3, 4)]
        for idx in indices_to_check:
            logger.debug(f"Checking index: {idx}")
            if 0 <= idx < len(wells):
                withdrawals = WellWithdrawal.objects.filter(target_well=wells[idx])
                if withdrawals:
                    return wells[idx], withdrawals
    return None, []


def get_new_plate_infos(experiment):
    plate_info = []
    plates = Plate.objects.filter(experiment=experiment)
    for plate in plates:
        logger.debug(f"Plate: {plate.barcode}")
        plate_details = PlateDetail.objects.get(pk=plate.id)
        measurement_labels = plate_details.measurement_labels
        measurement_timestamps = plate_details.measurement_timestamps

        if not measurement_labels or not measurement_timestamps:
            continue

        for label in measurement_labels:
            logger.info(f"Label: {label}")
            timestamps = measurement_timestamps.get(label, [])
            if not timestamps:
                continue
            for timestamp in timestamps:
                logger.info(f"Timestamp: {timestamp}")
                plate_info_obj = {
                    "measurement_label": label,
                    "measurement_timestamp": timestamp,
                    "replicate": "",
                    "cell_type": "",
                    "condition": "",
                }
                middle_well_index = int(
                    len(plate.wells.all()) // 2 + plate.dimension.cols // 2
                )  # we take one in the middle of the plate (the middle index plus the half number of columns), so we don't get a well filled from a control well.

                middle_well, withdrawals = find_withdrawal_well(
                    middle_well_index, plate.wells.all(), plate.dimension.cols
                )
                lib_plate = withdrawals[-1].well.plate if withdrawals else None
                plate_info_obj["plate_barcode"] = plate.barcode
                plate_info_obj["lib_plate_barcode"] = (
                    lib_plate.barcode if lib_plate else "NA"
                )
                plate_info.append(plate_info_obj)
    return plate_info


def prefill_plate_info(request):
    try:
        if request.method == "GET":
            experiment_id = request.GET.get("experiment_id")
            if not experiment_id:
                return JsonResponse({"error": "Experiment ID not provided"}, status=400)

            existing_plate_info = get_existing_plate_infos(experiment_id)
            if existing_plate_info:
                return JsonResponse({"plate_info": existing_plate_info}, status=200)

            experiment = get_object_or_404(Experiment, pk=experiment_id)
            new_plate_info = get_new_plate_infos(experiment)
            return JsonResponse({"plate_info": new_plate_info}, status=200)
    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
def save_plate_info(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            experiment_id = data.get("experiment_id")
            plate_info = data.get("plate_info")
            if not experiment_id:
                return JsonResponse({"error": "Experiment ID not provided"}, status=400)
            if not plate_info:
                return JsonResponse({"error": "Plate info not provided"}, status=400)

            experiment = Experiment.objects.get(pk=experiment_id)

            for item in plate_info:
                plate = Plate.objects.get(barcode=item["plate_barcode"])
                defaults = {
                    "lib_plate_barcode": item["lib_plate_barcode"],
                    "label": item["measurement_label"],
                    "replicate": item["replicate"],
                    "measurement_time": item["measurement_timestamp"],
                    "cell_type": item["cell_type"],
                    "condition": item["condition"],
                }
                PlateInfo.objects.update_or_create(
                    plate=plate, experiment=experiment, defaults=defaults
                )

            return JsonResponse({"status": "Plate info saved successfully"}, status=200)
    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)


def refresh(request):
    """
    Hard refresh on the ui side.
    When we delete measurements in admin, it will not be imeediately reflected on the UI,
    because we need to refresh materialized views.
    This function is triggerd by the refresh button on the ui.
    """
    try:
        if request.method == "GET":
            PlateDetail.refresh(concurrently=True)
            WellDetail.refresh(concurrently=True)
            ExperimentDetail.refresh(concurrently=True)
            return JsonResponse({"status": "Data refreshed successfully"}, status=200)
    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)
