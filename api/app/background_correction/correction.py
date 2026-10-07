"""
Corrects one measurement of a plate by its reference wells and saves the result
as a new measurement of the same plate, e.g. "Lum_CTG" -> "Lum_CTG_bc_N1_median".
Every time point is corrected by the reference wells of that time point.
"""

from rest_framework.exceptions import ValidationError

from background_correction.calculation import subtract_background
from background_correction.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
)
from core.models import Measurement, Plate


def corrected_label(label: str, reference_type: str, method: str) -> str:
    """("Lum_CTG", "N1", "median") -> "Lum_CTG_bc_N1_median" """
    return f"{label}_bc_{reference_type}_{method}"


def correct_plate(plate: Plate, label: str, reference_type: str, method: str) -> str:
    """Saves the corrected measurement and returns its label."""
    new_label = corrected_label(label, reference_type, method)
    check_label_length(new_label)
    values = values_by_time_point(plate, label)

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

    # A new correction with the same settings replaces the last one
    replace_measurements(plate, new_label, new_measurements)
    return new_label
