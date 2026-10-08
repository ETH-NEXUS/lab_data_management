"""
The activity of every well in percent between the controls of its plate, saved
as a new measurement, e.g. "Lum1" -> "Lum1_activity_N_P":

    %Activity = 100 * (value - median(P)) / (median(N) - median(P))

The negative control (N, nothing inhibits the kinase) is 100 %, the positive
control (P, fully inhibited) is 0 %. It is calculated from the raw readout,
separately for each time point. The negative control wells are left empty.
"""

import statistics

from rest_framework.exceptions import ValidationError

from core.models import Measurement, Plate
from plate_calculations.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
    values_of_wells,
    well_ids_of_type,
)


def activity_label(label: str, negative_type: str, positive_type: str) -> str:
    """("Lum1", "N", "P") -> "Lum1_activity_N_P" """
    return f"{label}_activity_{negative_type}_{positive_type}"


def percent_activity(
    values: dict[int, float],
    negative_well_ids: set[int],
    negative_median: float,
    positive_median: float,
) -> dict[int, float]:
    """
    The activity of every well but the negative controls, by well id.

    values = {11: 300.0, 21: 0.0, 31: 105.0}, negative = {11},
    negative_median = 200, positive_median = 10 -> {21: -5.26..., 31: 50.0}
    """
    activity = {}
    for well_id, value in values.items():
        if well_id not in negative_well_ids:
            activity[well_id] = (
                100 * (value - positive_median) / (negative_median - positive_median)
            )
    return activity


def activity_of_plate(
    plate: Plate, label: str, negative_type: str, positive_type: str
) -> str:
    """Saves the %Activity measurement and returns its label."""
    new_label = activity_label(label, negative_type, positive_type)
    check_label_length(new_label)
    values = values_by_time_point(plate, label)
    negative_well_ids = well_ids_of_type(plate, negative_type)
    positive_well_ids = well_ids_of_type(plate, positive_type)

    new_measurements = []
    for measured_at, well_values in values.items():
        negative_median = statistics.median(
            values_of_wells(
                plate, label, negative_type, negative_well_ids, well_values, measured_at
            )
        )
        positive_median = statistics.median(
            values_of_wells(
                plate, label, positive_type, positive_well_ids, well_values, measured_at
            )
        )
        # Without a difference between the controls there is no scale
        if negative_median == positive_median:
            raise ValidationError(
                f'The medians of the "{negative_type}" and "{positive_type}" wells are '
                f"the same at {measured_at}, so there is no %Activity to calculate."
            )
        activity = percent_activity(
            well_values, negative_well_ids, negative_median, positive_median
        )
        for well_id, value in activity.items():
            new_measurements.append(
                Measurement(
                    well_id=well_id,
                    label=new_label,
                    value=value,
                    measured_at=measured_at,
                )
            )

    # A new %Activity with the same controls replaces the last one
    replace_measurements(plate, new_label, new_measurements)
    return new_label
