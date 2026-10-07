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

from background_correction.plate_measurements import (
    check_label_length,
    replace_measurements,
    values_by_time_point,
)
from core.models import Measurement, Plate


def activity_label(label: str, negative_type: str, positive_type: str) -> str:
    """("Lum1", "N", "P") -> "Lum1_activity_N_P" """
    return f"{label}_activity_{negative_type}_{positive_type}"


def percent_activity(
    values: dict[int, float], negative_well_ids: set[int], positive_well_ids: set[int]
) -> dict[int, float]:
    """
    The activity of every well but the negative controls, by well id.

    values = {11: 100.0, 12: 300.0, 21: 0.0, 22: 20.0, 31: 55.0},
    negative = {11, 12} (median 200), positive = {21, 22} (median 10)
    -> {21: -5.26..., 22: 5.26..., 31: 23.68...}
    """
    negative = statistics.median(
        value for well_id, value in values.items() if well_id in negative_well_ids
    )
    positive = statistics.median(
        value for well_id, value in values.items() if well_id in positive_well_ids
    )
    activity = {}
    for well_id, value in values.items():
        if well_id not in negative_well_ids:
            activity[well_id] = 100 * (value - positive) / (negative - positive)
    return activity


def activity_of_plate(
    plate: Plate, label: str, negative_type: str, positive_type: str
) -> str:
    """Saves the %Activity measurement and returns its label."""
    new_label = activity_label(label, negative_type, positive_type)
    check_label_length(new_label)
    values = values_by_time_point(plate, label)

    well_ids = {}
    for well_type in [negative_type, positive_type]:
        well_ids[well_type] = set(
            plate.wells.filter(type__name=well_type).values_list("id", flat=True)
        )

    new_measurements = []
    for measured_at, well_values in values.items():
        for well_type in [negative_type, positive_type]:
            if not well_ids[well_type].intersection(well_values.keys()):
                raise ValidationError(
                    f'The plate {plate.barcode} has no "{well_type}" wells with a '
                    f'value of "{label}" measured at {measured_at}.'
                )
        negative_values = [
            well_values[i] for i in well_ids[negative_type] if i in well_values
        ]
        positive_values = [
            well_values[i] for i in well_ids[positive_type] if i in well_values
        ]
        # Without a difference between the controls there is no scale
        if statistics.median(negative_values) == statistics.median(positive_values):
            raise ValidationError(
                f'The medians of the "{negative_type}" and "{positive_type}" wells are '
                f"the same at {measured_at}, so there is no %Activity to calculate."
            )
        activity = percent_activity(
            well_values, well_ids[negative_type], well_ids[positive_type]
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
