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

from analysis.controls import check_chosen_controls, control_warnings, rename_controls
from core.models import Experiment, Measurement
from core.views.plate_info import get_existing_plate_infos
from ldm.ldm import get_experiment_measurements

# The columns of get_experiment_measurements(type="main") in ldm/ldm.py. The
# "chemical" export has the same rows with the compound data added, so main_info
# is taken from it instead of reading every well a second time (about 20 s for a
# screen of 26 plates). A test checks that both give the same main_info.
MAIN_COLUMNS = [
    "unique_identifier",
    "well_coordinate",
    "value",
    "plate",
    "plate_row",
    "plate_column",
    "control",
    "measurement",
    "is_invalid",
    "measured_at",
    "compound",
]

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
    experiment: Experiment,
    label: str,
    folder: str,
    conditions: list[str],
    controls: dict,
) -> tuple[dict, list[str]]:
    """
    Writes the three files into `folder`. `conditions` are the ones the report
    compares (none for a single analysis, two for a selectivity analysis).
    `controls` are the well types of the controls, e.g. {"positive": "P1",
    "negative": "N1"}; they are written as "P" and "N", the names R knows.
    Returns the paths of the files, by the name of the report parameter that
    reads them, and warnings about the data:
    ({"path_data": "/x/main_info.csv", "path_lib": "/x/chemical_info.csv",
      "path_meta": "/x/experiment_data.csv"},
     ['Plate SP_3 has no positive control wells ("P1"), so ...'])

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
        text = (
            f'The experiment "{experiment.name}" has no plate information for the '
            f'measurement "{label}". The report needs it for every plate (library '
            "plate, replicate, cell type, condition). Add it on the experiment page "
            'with "add experiment data" and save it.'
        )
        if labels:
            text += f" Plate information exists for: {quoted(labels)}."
        else:
            text += " The experiment has no plate information for any measurement yet."
        raise CommandError(text)
    check_conditions(conditions, plate_infos)

    chemical_info = get_experiment_measurements(
        experiment.name, label, "chemical", csv=True, experiment_id=experiment.id
    )
    if chemical_info.empty:
        labels = (
            Measurement.objects.filter(well__plate__experiment=experiment)
            .values_list("label", flat=True)
            .distinct()
        )
        raise CommandError(
            f'The experiment "{experiment.name}" has no measurements "{label}". '
            f"Its measurements are: {quoted(sorted(labels)) or 'none'}."
        )
    chemical_info, repeated_warnings = remove_repeated_readings(chemical_info, label)
    positive, negative = controls["positive"], controls["negative"]
    check_chosen_controls(chemical_info, positive, negative)
    chemical_info = rename_controls(chemical_info, positive, negative)
    main_info = chemical_info[MAIN_COLUMNS]
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
    warnings = repeated_warnings + plate_warnings(main_info, experiment_data)
    warnings += control_warnings(main_info, positive, negative)
    return paths, warnings


def remove_repeated_readings(
    chemical_info: pd.DataFrame, label: str
) -> tuple[pd.DataFrame, list[str]]:
    """
    The report needs one reading per well. A measurement file that was mapped
    twice gives each well a second row with the same value and another time
    (seen on prod: experiments 83, 86 and 110); the copy is removed, with a
    warning. Two readings with different values cannot be told apart by the
    report, so they stop the analysis.

    Returned data example:
    (the rows without the copies, ['The measurement "Lum" was mapped twice: ...'])
    """
    columns_without_time = [
        column for column in chemical_info.columns if column != "measured_at"
    ]
    without_copies = chemical_info.drop_duplicates(subset=columns_without_time)
    copies = len(chemical_info) - len(without_copies)

    repeated = without_copies[without_copies["unique_identifier"].duplicated()]
    if not repeated.empty:
        example = repeated.iloc[0]
        raise CommandError(
            f'The measurement "{label}" has more than one reading with different values '
            f"on {repeated['unique_identifier'].nunique()} well(s) (e.g. well "
            f"{example['well_coordinate']} of plate {example['plate']}). The report can "
            "use one reading per well only: map the readings with different measurement "
            "names, or keep only the reading that belongs to the analysis."
        )

    warnings = []
    if copies > 0:
        warnings.append(
            f'The measurement "{label}" was mapped more than once: {copies} readings '
            "are copies of another reading of the same well with the same value. "
            "Each well is used once in the report."
        )
    return without_copies, warnings


def quoted(names: list[str]) -> str:
    """
    Names in quotes, so a space at the end is seen.
    ["Lum", "Log "] -> '"Lum", "Log "'
    """
    return ", ".join(f'"{name}"' for name in names)


def plate_warnings(main_info: pd.DataFrame, experiment_data: pd.DataFrame) -> list[str]:
    """Plates with measurements but no plate information, and the other way round."""
    measured = set(main_info["plate"])
    described = set(experiment_data["plate"])
    warnings = []
    without_info = sorted(measured - described)
    if without_info:
        warnings.append(
            "These plates have measurements but no plate information, so the report "
            "has no library plate, replicate, cell type or condition for them: "
            f"{', '.join(without_info)}. Add it with \"add experiment data\"."
        )
    without_measurements = sorted(described - measured)
    if without_measurements:
        warnings.append(
            "These plates have plate information but no measurements of this label, "
            f"so they are not in the report: {', '.join(without_measurements)}."
        )
    return warnings


def check_conditions(conditions: list[str], plate_infos: list[dict]) -> None:
    """
    A selectivity report compares two conditions of the plate information. With a
    condition that is not there, R fails with an error that does not say why.
    """
    known = sorted({str(plate_info["condition"]) for plate_info in plate_infos})
    unknown = [condition for condition in conditions if condition not in known]
    if unknown and known == [""]:
        raise CommandError(
            "A selectivity analysis compares two conditions of the plate information, "
            "but no plate of this measurement has a condition yet. Fill in the "
            'column "Condition" with "add experiment data" and save it.'
        )
    if unknown:
        raise CommandError(
            "A selectivity analysis compares two conditions of the plate information, "
            f"but these are not in it: {quoted(unknown)}. The conditions of this "
            f"measurement are: {quoted(known)}. Type two of them exactly like this, "
            'or correct the conditions with "add experiment data".'
        )
