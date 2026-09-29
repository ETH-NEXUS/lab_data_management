"""
Plate information of an experiment: prefilling and saving it.
"""

from django.shortcuts import get_object_or_404
from django.http import JsonResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated

from helpers.logger import logger
from ..utils.plates.plate_infos import with_saved_values
from ..models import (
    Plate,
    Experiment,
    PlateDetail,
    PlateInfo,
)


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


def find_well_with_donors(start_index, wells, total_columns):
    """
    If the middle well is empty, we look for the closest well with withdrawals up and down.
    """
    if len(wells[start_index].donors.all()) > 0:
        return wells[start_index]
    indices_to_check = [start_index + i * total_columns for i in range(-4, 5)]
    for idx in indices_to_check:
        if idx < 0 or idx >= len(wells):
            continue
        if len(wells[idx].donors.all()) > 0:
            return wells[idx]

    return None


def get_new_plate_infos(experiment):
    plate_info = []
    plates = Plate.objects.filter(experiment=experiment)
    for plate in plates:
        logger.info(f"Plate: {plate.barcode}")
        # A plate mapped after the last refresh of the materialized views is not there yet
        plate_details = PlateDetail.objects.filter(pk=plate.id).first()
        if plate_details is None:
            continue
        measurement_labels = plate_details.measurement_labels
        measurement_timestamps = plate_details.measurement_timestamps

        if not measurement_labels or not measurement_timestamps:
            continue

        for label in measurement_labels:
            logger.info(f"Label: {label}")
            timestamps = measurement_timestamps.get(label, [])
            if not timestamps:
                continue
            # One row per plate and label, because the plate information is saved
            # per plate and label (save_plate_info). A label that was read several
            # times gets the time of its first reading, e.g. of
            # ["2025-05-16T10:00:00+02:00", "2025-05-17T10:00:00+02:00"] the first.
            timestamp = min(timestamps)
            logger.info(f"Timestamp: {timestamp}")
            plate_info_obj = {
                "measurement_label": label,
                "measurement_timestamp": timestamp,
                "replicate": "",
                "cell_type": "",
                "condition": "",
            }
            well_to_check_idx = int(
                len(plate.wells.all()) // 2 + plate.dimension.cols // 2
            )  # we take one in the middle of the plate (the middle index plus the half number of columns), so we don't get a well filled from a control well.

            well_with_donors = find_well_with_donors(
                well_to_check_idx,
                plate.wells.all().order_by("position"),
                plate.dimension.cols,
            )

            lib_plate = None
            if well_with_donors:
                lib_plate = well_with_donors.donors.all().first().well.plate
            plate_info_obj["plate_barcode"] = plate.barcode
            plate_info_obj["lib_plate_barcode"] = (
                lib_plate.barcode if lib_plate else "NA"
            )
            plate_info.append(plate_info_obj)

    return plate_info


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def prefill_plate_info(request):
    """
    The plate information of an experiment for the form: one row per plate and
    measurement label, with the values that were saved before.
    """
    experiment_id = request.GET.get("experiment_id")
    if not experiment_id:
        return JsonResponse({"error": "Experiment ID not provided"}, status=400)

    experiment = get_object_or_404(Experiment, pk=experiment_id)
    plate_info = with_saved_values(
        get_new_plate_infos(experiment), get_existing_plate_infos(experiment_id)
    )
    return JsonResponse({"plate_info": plate_info}, status=200)


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def save_plate_info(request):
    """Stores the plate information of an experiment."""
    experiment_id = request.data.get("experiment_id")
    plate_info = request.data.get("plate_info")
    if not experiment_id:
        return JsonResponse({"error": "Experiment ID not provided"}, status=400)
    if not plate_info:
        return JsonResponse({"error": "Plate info not provided"}, status=400)

    experiment = Experiment.objects.get(pk=experiment_id)
    for item in plate_info:
        plate = Plate.objects.get(barcode=item["plate_barcode"])
        # One row per plate and measurement label: a plate measured twice (e.g. Lum
        # and Fluo) has two rows, and the second must not overwrite the first
        PlateInfo.objects.update_or_create(
            plate=plate,
            experiment=experiment,
            label=item["measurement_label"],
            defaults={
                "lib_plate_barcode": item["lib_plate_barcode"],
                "replicate": item["replicate"],
                "measurement_time": item["measurement_timestamp"],
                "cell_type": item["cell_type"],
                "condition": item["condition"],
            },
        )
    return JsonResponse({"status": "Plate info saved successfully"}, status=200)
