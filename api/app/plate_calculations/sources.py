"""
Where the measurements of a plate come from: the file they were imported from,
or none (calculated on the plate page or with the measurement calculator of the
experiment). The plate page names the file below the heatmap.
"""

import os

from core.models import Measurement, Plate


def measurement_sources(plate: Plate) -> dict[str, str | None]:
    """
    The file name of every imported measurement of the plate, None for the others.

    {"Lum1": "093026-154654_RKS_300926_1.asc", "Lum1_log10": None}
    """
    sources: dict[str, str | None] = {}
    rows = (
        Measurement.objects.filter(well__plate=plate)
        .values_list("label", "measurement_assignment__filename")
        .distinct()
    )
    for label, file_path in rows:
        if file_path:
            # The full path of the server says nothing to the lab, the name does
            sources[label] = os.path.basename(file_path)
        elif label not in sources:
            sources[label] = None
    return sources
