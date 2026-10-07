"""
Reading and saving the measurements of one plate, shared by the calculations of
this app (background correction, log10).
"""

from django.db import transaction
from rest_framework.exceptions import ValidationError

from core.models import (
    ExperimentDetail,
    Measurement,
    Plate,
    PlateDetail,
    WellDetail,
)

LABEL_MAX_LENGTH = Measurement._meta.get_field("label").max_length


def check_label_length(new_label: str) -> None:
    if len(new_label) > LABEL_MAX_LENGTH:
        raise ValidationError(
            f'The name of the new measurement "{new_label}" is longer than '
            f"{LABEL_MAX_LENGTH} characters."
        )


def values_by_time_point(plate: Plate, label: str) -> dict:
    """
    The values of one measurement of the plate, per time point and well id.
    A plate without the measurement is refused.

    {datetime(2025, 5, 16, 10, 0): {11: 10.0, 12: 4.0}}
    """
    measurements = Measurement.objects.filter(well__plate=plate, label=label)
    values: dict = {}
    for measurement in measurements:
        if measurement.measured_at not in values:
            values[measurement.measured_at] = {}
        values[measurement.measured_at][measurement.well_id] = measurement.value
    if not values:
        raise ValidationError(
            f'The plate {plate.barcode} has no measurement "{label}".'
        )
    return values


def replace_measurements(
    plate: Plate, new_label: str, new_measurements: list[Measurement]
) -> None:
    """
    Saves the new measurement of the plate instead of an earlier one of the same
    label, and refreshes the views the plate page reads.
    """
    with transaction.atomic():
        Measurement.objects.filter(well__plate=plate, label=new_label).delete()
        Measurement.objects.bulk_create(new_measurements)

    PlateDetail.refresh(concurrently=True)
    WellDetail.refresh(concurrently=True)
    ExperimentDetail.refresh(concurrently=True)
