"""
The three input files of the statistics reports (api/app/statistics/*.qmd) for one
experiment and one measurement label. They have the columns the statistics group
used to download from LDM by hand, and that SLmisc.R reads:

- main_info.csv        one row per well: unique_identifier, plate, plate_row,
                       plate_column, control, value, ...
- chemical_info.csv    the same rows, plus the compound data of the library
- experiment_data.csv  one row per plate: measurement_label, measurement_timestamp,
                       plate, lib_plate_barcode, replicate, cell_type, condition
"""

import os

import pandas as pd
from django.core.management.base import CommandError

from core.models import Experiment
from core.views.plate_info import get_existing_plate_infos
from ldm.ldm import get_experiment_measurements

# Same columns and order as downloadCsv() in ui/app/pages/add_data/[experiment_id].vue
EXPERIMENT_DATA_COLUMNS = [
    "measurement_label",
    "measurement_timestamp",
    "plate",
    "lib_plate_barcode",
    "replicate",
    "cell_type",
    "condition",
]


def write_input_files(
    experiment: Experiment, label: str, folder: str, conditions: list[str]
) -> tuple[dict, list[str]]:
    """
    Writes the three files into `folder`. `conditions` are the ones the report
    compares (none for a single analysis, two for a selectivity analysis).
    Returns the paths of the files, by the name of the report parameter that
    reads them, and warnings about the data:
    ({"path_data": "/x/main_info.csv", "path_lib": "/x/chemical_info.csv",
      "path_meta": "/x/experiment_data.csv"},
     ["Plate SP_3 has no P control wells: ..."])

    A report reads one measurement label only: it joins the wells and the plate
    information by plate, so a second label would pair each value with both labels.
    """
    all_plate_infos = get_existing_plate_infos(experiment.id)
    plate_infos = [
        plate_info
        for plate_info in all_plate_infos
        if plate_info["measurement_label"] == label
    ]
    if not plate_infos:
        labels = sorted({info["measurement_label"] for info in all_plate_infos})
        raise CommandError(
            f'The experiment "{experiment.name}" has no plate information for the '
            f'measurement "{label}". Add it with "add experiment data". '
            f"Plate information exists for: {', '.join(labels) or 'no measurement'}."
        )
    check_conditions(conditions, plate_infos)

    main_info = get_experiment_measurements(experiment.name, label, "main", csv=True)
    if main_info.empty:
        raise CommandError(
            f'The experiment "{experiment.name}" has no measurements "{label}".'
        )
    chemical_info = get_experiment_measurements(
        experiment.name, label, "chemical", csv=True
    )
    experiment_data = pd.DataFrame(plate_infos).rename(
        columns={"plate_barcode": "plate"}
    )

    paths = {
        "path_data": os.path.join(folder, "main_info.csv"),
        "path_lib": os.path.join(folder, "chemical_info.csv"),
        "path_meta": os.path.join(folder, "experiment_data.csv"),
    }
    main_info.to_csv(paths["path_data"], index=False)
    chemical_info.to_csv(paths["path_lib"], index=False)
    experiment_data[EXPERIMENT_DATA_COLUMNS].to_csv(paths["path_meta"], index=False)
    return paths, control_warnings(main_info)


def check_conditions(conditions: list[str], plate_infos: list[dict]) -> None:
    """
    A selectivity report compares two conditions of the plate information. With a
    condition that is not there, R fails with an error that does not say why.
    """
    known = sorted({str(plate_info["condition"]) for plate_info in plate_infos})
    unknown = [condition for condition in conditions if condition not in known]
    if unknown:
        raise CommandError(
            f"Condition not in the plate information: {', '.join(unknown)}. "
            f"The conditions of this measurement are: {', '.join(known)}."
        )


def control_warnings(main_info: pd.DataFrame) -> list[str]:
    """
    The report normalizes every plate by its negative (N) and positive (P) control
    wells; a plate without them is left out of the results (seen on experiment 105:
    plate 20250513SP_29 has no P wells and no results).
    """
    warnings = []
    for plate, controls in main_info.groupby("plate")["control"]:
        missing = [control for control in ["N", "P"] if control not in set(controls)]
        if missing:
            warnings.append(
                f"Plate {plate} has no {' and no '.join(missing)} control wells, so "
                "the report cannot normalize it and leaves it out of the results."
            )
    return warnings
