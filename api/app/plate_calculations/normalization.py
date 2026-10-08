"""
Normalization of one measurement of a plate between its controls, as in the R
report of the lab (Michael's analysis), saved as two new measurements:

    %Inhibition = (log10(x) - median(log10 N)) / (median(log10 P) - median(log10 N))
    %Activity   = (1 - %Inhibition) * 100

The negative control (N, nothing inhibits the enzyme) has an inhibition of 0 and
an activity of 100, the positive control (P, fully inhibited) an inhibition of 1
and an activity of 0. The %Inhibition is a fraction like in the report, the
%Activity is in percent, as the lab asked for. Every time point is
normalized by the controls of that time point. A well with a value of 0 or below
has no log10 and is left empty in every read, as in the log10 calculation.
"""

import math
import statistics

from rest_framework.exceptions import ValidationError

from core.models import Measurement, Plate
from plate_calculations.log_transform import wells_without_log10
from plate_calculations.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
    values_of_wells,
    well_ids_of_type,
)


def inhibition_label(label: str, negative_type: str, positive_type: str) -> str:
    """("Lum1", "N", "P") -> "Lum1_inhibition_N_P" """
    return f"{label}_inhibition_{negative_type}_{positive_type}"


def activity_label(label: str, negative_type: str, positive_type: str) -> str:
    """("Lum1", "N", "P") -> "Lum1_activity_N_P" """
    return f"{label}_activity_{negative_type}_{positive_type}"


def inhibition(
    logs: dict[int, float], negative_median: float, positive_median: float
) -> dict[int, float]:
    """
    The inhibition of every well from the log10 of its value, by well id.

    logs = {11: 3.0, 21: 2.0, 31: 2.5}, negative_median = 3.0,
    positive_median = 2.0 -> {11: 0.0, 21: 1.0, 31: 0.5}
    """
    result = {}
    for well_id, log in logs.items():
        result[well_id] = (log - negative_median) / (positive_median - negative_median)
    return result


def normalize_plate(
    plate: Plate, label: str, negative_type: str, positive_type: str
) -> tuple[str, str, int]:
    """
    Saves the %Inhibition and %Activity measurements and returns their labels and
    how many wells were left empty, e.g. ("Lum1_inhibition_N_P", "Lum1_activity_N_P", 0).
    """
    new_inhibition_label = inhibition_label(label, negative_type, positive_type)
    new_activity_label = activity_label(label, negative_type, positive_type)
    check_label_length(new_inhibition_label)
    values = values_by_time_point(plate, label)
    skipped_wells = wells_without_log10(values)
    negative_well_ids = well_ids_of_type(plate, negative_type)
    positive_well_ids = well_ids_of_type(plate, positive_type)

    inhibition_measurements = []
    activity_measurements = []
    for measured_at, well_values in values.items():
        logs = {
            well_id: math.log10(value)
            for well_id, value in well_values.items()
            if well_id not in skipped_wells
        }
        negative_median = statistics.median(
            values_of_wells(
                plate, label, negative_type, negative_well_ids, logs, measured_at
            )
        )
        positive_median = statistics.median(
            values_of_wells(
                plate, label, positive_type, positive_well_ids, logs, measured_at
            )
        )
        # Without a difference between the controls there is no scale
        if negative_median == positive_median:
            raise ValidationError(
                f'The medians of the "{negative_type}" and "{positive_type}" wells are '
                f"the same at {measured_at}, so the values cannot be normalized."
            )
        for well_id, value in inhibition(
            logs, negative_median, positive_median
        ).items():
            inhibition_measurements.append(
                Measurement(
                    well_id=well_id,
                    label=new_inhibition_label,
                    value=value,
                    measured_at=measured_at,
                )
            )
            activity_measurements.append(
                Measurement(
                    well_id=well_id,
                    label=new_activity_label,
                    value=(1 - value) * 100,
                    measured_at=measured_at,
                )
            )

    # A new normalization with the same controls replaces the last one
    replace_measurements(
        plate,
        {
            new_inhibition_label: inhibition_measurements,
            new_activity_label: activity_measurements,
        },
    )
    return new_inhibition_label, new_activity_label, len(skipped_wells)
