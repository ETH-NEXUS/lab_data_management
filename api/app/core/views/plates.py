"""
Plates, including their measurements and time points.
"""

from datetime import datetime
import traceback
from django.db.models import Prefetch, Q
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.exceptions import ValidationError
from rest_framework.response import Response
from compoundlib.serializers import SimpleCompoundLibrarySerializer
from helpers.logger import logger
from ..models import (
    Well,
    Plate,
    Measurement,
    WellWithdrawal,
    WellCompound,
    Experiment,
    MeasurementFeature,
    PlateDetail,
    WellDetail,
    ExperimentDetail,
)
from ..serializers import (
    PlateSerializer,
    SimpleExperimentSerializer,
    SimplePlateTemplateSerializer,
)
from .plate_archive import PlateArchiveMixin
from ..utils.plates.archive_guard import (
    ensure_plate_can_be_changed,
    is_archived_library_plate,
)


GLOBAL_NOW = datetime.now().replace(microsecond=0)


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


class PlateViewSet(PlateArchiveMixin, viewsets.ModelViewSet):
    def get_serializer_class(self):
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

    # Archived plates cannot be changed or deleted through the API
    # (see core/utils/plates/archive_guard.py). Archiving itself is a separate action.
    def perform_update(self, serializer):
        ensure_plate_can_be_changed(serializer.instance)
        super().perform_update(serializer)

    def perform_destroy(self, instance):
        ensure_plate_can_be_changed(instance)
        super().perform_destroy(instance)

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

    @action(detail=True, methods=["post"])
    def apply_template(self, request, pk=None):
        """Applies a template plate"""
        apply_to_all_experiment_plates = request.data.get(
            "apply_to_all_experiment_plates"
        )
        template_plate_id = request.data.get("template")
        if template_plate_id is None:
            raise ValidationError(
                {"template": ["This field is required."]}, code="invalid"
            )

        template_plate = Plate.objects.get(pk=template_plate_id)
        plate = self.get_object()
        ensure_plate_can_be_changed(plate)

        if apply_to_all_experiment_plates:
            plates = Plate.objects.filter(experiment=plate.experiment)
            for _plate in plates:
                # Archived library plates are left as they are, instead of making the
                # whole request fail because of a plate nobody asked to change.
                if is_archived_library_plate(_plate):
                    continue
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
        result = None
        try:

            import math  # don't remove this import!!!!!

            print(
                math
            )  # don't remove this print!!!!! we need it so that the IDE don't remove the import by formating code
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
        except Exception as e:
            traceback.print_exc()
            logger.critical(
                f"Error evaluating expression: {e}. It can be that the value in the expression {new_expression} is <= 0. Setting result to 0"
            )
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
