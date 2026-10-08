"""
The log10 of one measurement of a plate, saved as a new measurement of the same
plate, e.g. "Lum1" -> "Lum1_log10". Readouts spread over several orders of
magnitude, and on the log scale the heatmap shows the differences between them.
A value of 0 or below has no log10: its well is left empty.
"""

import math

from rest_framework.exceptions import ValidationError

from core.models import Measurement, Plate
from plate_calculations.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
)


def log10_label(label: str) -> str:
    """ "Lum1" -> "Lum1_log10" """
    return f"{label}_log10"


def wells_without_log10(values: dict) -> set[int]:
    """
    The wells with a value of 0 or below in any read. They are left empty in
    every read: the plate page lists the values of a well by the order of the
    reads, so a gap in one read would move the later values to the wrong read.

    {t1: {11: 1000.0, 12: 0.0}, t2: {11: 10.0, 12: 5.0}} -> {12}
    """
    wells = set()
    for well_values in values.values():
        for well_id, value in well_values.items():
            if value <= 0:
                wells.add(well_id)
    return wells


def log10_of_plate(plate: Plate, label: str) -> tuple[str, int]:
    """
    Saves the log10 measurement and returns its label and how many wells were
    left empty because a value is 0 or below, e.g. ("Lum1_log10", 2).
    """
    new_label = log10_label(label)
    check_label_length(new_label)
    values = values_by_time_point(plate, label)
    skipped_wells = wells_without_log10(values)

    new_measurements = []
    for measured_at, well_values in values.items():
        for well_id, value in well_values.items():
            if well_id in skipped_wells:
                continue
            new_measurements.append(
                Measurement(
                    well_id=well_id,
                    label=new_label,
                    value=math.log10(value),
                    measured_at=measured_at,
                )
            )

    if not new_measurements:
        raise ValidationError(
            f'Every well of "{label}" on the plate {plate.barcode} has a value of 0 '
            "or below, so there is no log10 to show."
        )

    # A new log10 of the same measurement replaces the last one
    replace_measurements(plate, {new_label: new_measurements})
    return new_label, len(skipped_wells)
