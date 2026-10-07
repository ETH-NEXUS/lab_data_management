"""
The background correction of one plate at one time point: the median (or mean)
of the reference wells is subtracted from every other well. The reference wells
are left out of the result, they are not needed after the correction.
"""

import statistics

# The ways to get the background from the values of the reference wells
METHODS = {
    "median": statistics.median,
    "mean": statistics.mean,
}


def subtract_background(
    values: dict[int, float], reference_well_ids: set[int], method: str
) -> dict[int, float]:
    """
    The corrected values of all wells that are not reference wells.

    values = {11: 10.0, 12: 4.0, 13: 2.0, 14: 6.0}, reference_well_ids = {13, 14},
    method = "median" -> background 4.0 -> {11: 6.0, 12: 0.0}
    """
    reference_values = []
    for well_id, value in values.items():
        if well_id in reference_well_ids:
            reference_values.append(value)

    background = METHODS[method](reference_values)

    corrected_values = {}
    for well_id, value in values.items():
        if well_id not in reference_well_ids:
            corrected_values[well_id] = value - background
    return corrected_values
