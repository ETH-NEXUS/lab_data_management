"""
Background correction of one measurement of a plate: the median (or mean) of the
reference wells is subtracted from every other well, separately for each time
point, and saved as a new measurement, e.g. "Lum_CTG" -> "Lum_CTG_bc_N1_median".
The reference wells are left empty, they are not needed after the correction.
"""

import statistics

from core.models import Measurement, Plate
from plate_calculations.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
    values_of_wells,
    well_ids_of_type,
)

# The ways to get the background from the values of the reference wells
METHODS = {
    "median": statistics.median,
    "mean": statistics.mean,
}


def corrected_label(label: str, reference_type: str, method: str) -> str:
    """("Lum_CTG", "N1", "median") -> "Lum_CTG_bc_N1_median" """
    return f"{label}_bc_{reference_type}_{method}"


def subtract_background(
    values: dict[int, float], reference_well_ids: set[int], background: float
) -> dict[int, float]:
    """
    The corrected values of all wells that are not reference wells.

    values = {11: 10.0, 12: 4.0, 13: 2.0}, reference_well_ids = {13},
    background = 4.0 -> {11: 6.0, 12: 0.0}
    """
    corrected_values = {}
    for well_id, value in values.items():
        if well_id not in reference_well_ids:
            corrected_values[well_id] = value - background
    return corrected_values


def correct_plate(plate: Plate, label: str, reference_type: str, method: str) -> str:
    """Saves the corrected measurement and returns its label."""
    new_label = corrected_label(label, reference_type, method)
    check_label_length(new_label)
    values = values_by_time_point(plate, label)
    reference_well_ids = well_ids_of_type(plate, reference_type)

    new_measurements = []
    for measured_at, well_values in values.items():
        reference_values = values_of_wells(
            plate, label, reference_type, reference_well_ids, well_values, measured_at
        )
        background = METHODS[method](reference_values)
        corrected_values = subtract_background(
            well_values, reference_well_ids, background
        )
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
