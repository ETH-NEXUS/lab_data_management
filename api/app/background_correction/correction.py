"""
Corrects one measurement of a plate by its reference wells and saves the result
as a new measurement of the same plate, e.g. "Lum_CTG" -> "Lum_CTG_bc_N1_median".
Every time point is corrected by the reference wells of that time point.
"""

from django.db import transaction
from rest_framework.exceptions import ValidationError

from background_correction.calculation import subtract_background
from core.models import (
    ExperimentDetail,
    Measurement,
    Plate,
    PlateDetail,
    WellDetail,
)

LABEL_MAX_LENGTH = Measurement._meta.get_field("label").max_length


def corrected_label(label: str, reference_type: str, method: str) -> str:
    """("Lum_CTG", "N1", "median") -> "Lum_CTG_bc_N1_median" """
    return f"{label}_bc_{reference_type}_{method}"


def values_by_time_point(plate: Plate, label: str) -> dict:
    """
    The values of one measurement of the plate, per time point and well id.

    {datetime(2025, 5, 16, 10, 0): {11: 10.0, 12: 4.0}}
    """
    measurements = Measurement.objects.filter(well__plate=plate, label=label)
    values = {}
    for measurement in measurements:
        if measurement.measured_at not in values:
            values[measurement.measured_at] = {}
        values[measurement.measured_at][measurement.well_id] = measurement.value
    return values


def correct_plate(plate: Plate, label: str, reference_type: str, method: str) -> str:
    """Saves the corrected measurement and returns its label."""
    new_label = corrected_label(label, reference_type, method)
    if len(new_label) > LABEL_MAX_LENGTH:
        raise ValidationError(
            f'The name of the corrected measurement "{new_label}" is longer than '
            f"{LABEL_MAX_LENGTH} characters."
        )

    values = values_by_time_point(plate, label)
    if not values:
        raise ValidationError(
            f'The plate {plate.barcode} has no measurement "{label}".'
        )

    reference_well_ids = set(
        plate.wells.filter(type__name=reference_type).values_list("id", flat=True)
    )

    new_measurements = []
    for measured_at, well_values in values.items():
        # Without reference values there is no background to subtract
        if not reference_well_ids.intersection(well_values.keys()):
            raise ValidationError(
                f'The plate {plate.barcode} has no "{reference_type}" wells with a '
                f'value of "{label}" measured at {measured_at}.'
            )
        corrected_values = subtract_background(well_values, reference_well_ids, method)
        for well_id, value in corrected_values.items():
            new_measurements.append(
                Measurement(
                    well_id=well_id,
                    label=new_label,
                    value=value,
                    measured_at=measured_at,
                )
            )

    with transaction.atomic():
        # A new correction with the same settings replaces the last one
        Measurement.objects.filter(well__plate=plate, label=new_label).delete()
        Measurement.objects.bulk_create(new_measurements)

    PlateDetail.refresh(concurrently=True)
    WellDetail.refresh(concurrently=True)
    ExperimentDetail.refresh(concurrently=True)
    return new_label
