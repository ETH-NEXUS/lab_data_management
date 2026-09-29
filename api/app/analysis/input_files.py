"""
The three input files of the statistics reports (api/app/statistics/*.qmd) for one
experiment and one measurement label. They have the columns the statistics group
used to download from LDM by hand, and that SLmisc.R reads:

- main_info.csv        one row per well: unique_identifier, plate, plate_row,
                       plate_column, control, value, ...
- chemical_info.csv    the same rows, plus the compound data of the library
- experiment_data.csv  one row per plate: measurement_label, measurement_timestamp,
                       plate, lib_plate_barcode, replicate, cell_type, condition

The data is checked first, so a problem is explained here and not by an R error
in the middle of the report.
"""

import os

import pandas as pd
from django.core.management.base import CommandError

from core.models import Experiment, Measurement
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

# The report normalizes every plate by its negative and positive control wells
CONTROL_TYPES = {"N": "negative", "P": "positive"}


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
     ["Plate SP_3 has no positive (P) control wells: ..."])

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
        raise CommandError(missing_plate_info_text(experiment, label, all_plate_infos))
    check_conditions(conditions, plate_infos)

    main_info = get_experiment_measurements(experiment.name, label, "main", csv=True)
    if main_info.empty:
        raise CommandError(missing_measurements_text(experiment, label))
    experiment_data = pd.DataFrame(plate_infos).rename(
        columns={"plate_barcode": "plate"}
    )
    warnings = check_controls(main_info) + check_plates(main_info, experiment_data)

    chemical_info = get_experiment_measurements(
        experiment.name, label, "chemical", csv=True
    )

    paths = {
        "path_data": os.path.join(folder, "main_info.csv"),
        "path_lib": os.path.join(folder, "chemical_info.csv"),
        "path_meta": os.path.join(folder, "experiment_data.csv"),
    }
    main_info.to_csv(paths["path_data"], index=False)
    chemical_info.to_csv(paths["path_lib"], index=False)
    experiment_data[EXPERIMENT_DATA_COLUMNS].to_csv(paths["path_meta"], index=False)
    return paths, warnings


def missing_plate_info_text(
    experiment: Experiment, label: str, all_plate_infos: list[dict]
) -> str:
    """Why the report cannot start without plate information, and what to do."""
    text = (
        f'The experiment "{experiment.name}" has no plate information for the '
        f'measurement "{label}". The report needs it for every plate (library '
        "plate, replicate, cell type, condition). Add it on the experiment page "
        'with "add experiment data".'
    )
    other_labels = sorted({info["measurement_label"] for info in all_plate_infos})
    if other_labels:
        text += f" Plate information exists for: {', '.join(other_labels)}."
    return text


def missing_measurements_text(experiment: Experiment, label: str) -> str:
    """Which measurement labels the experiment does have."""
    labels = (
        Measurement.objects.filter(well__plate__experiment=experiment)
        .values_list("label", flat=True)
        .distinct()
    )
    text = f'The experiment "{experiment.name}" has no measurements "{label}".'
    if labels:
        text += f" Its measurements are: {', '.join(sorted(labels))}."
    else:
        text += " It has no measurements at all yet."
    return text


def check_conditions(conditions: list[str], plate_infos: list[dict]) -> None:
    """
    A selectivity report compares two conditions of the plate information. With a
    condition that is not there, R fails with an error that does not say why.
    """
    known = sorted({str(plate_info["condition"]) for plate_info in plate_infos})
    unknown = [condition for condition in conditions if condition not in known]
    if unknown:
        raise CommandError(
            "Not in the plate information of this measurement: "
            f"{', '.join(repr(c) for c in unknown)}. Its conditions are: "
            f"{', '.join(repr(c) for c in known)}. Choose two of them, or correct "
            'the conditions with "add experiment data".'
        )


def check_controls(main_info: pd.DataFrame) -> list[str]:
    """
    The report normalizes every plate by its N and P control wells. Without them
    in any plate it cannot run at all; a single plate without them is left out
    (seen on experiment 105: plate 20250513SP_29 has no P wells and no results).
    """
    controls_per_plate = main_info.groupby("plate")["control"].apply(set)
    all_controls = set().union(*controls_per_plate)
    for control, name in CONTROL_TYPES.items():
        if control not in all_controls:
            raise CommandError(
                f"No plate of this measurement has {name} ({control}) control wells, "
                "and the report normalizes every plate by the N and P controls. "
                f"The well types found are: {', '.join(sorted(all_controls))}. "
                "Check the well types of the plate layout."
            )

    warnings = []
    for plate, controls in controls_per_plate.items():
        missing = [
            f"{CONTROL_TYPES[control]} ({control})"
            for control in CONTROL_TYPES
            if control not in controls
        ]
        if missing:
            warnings.append(
                f"Plate {plate} has no {' and no '.join(missing)} control wells: "
                "the report cannot normalize it, so its wells are not in the results."
            )
    return warnings


def check_plates(main_info: pd.DataFrame, experiment_data: pd.DataFrame) -> list[str]:
    """Plates with measurements but without plate information, and the other way round."""
    measured = set(main_info["plate"])
    described = set(experiment_data["plate"])
    warnings = []
    without_info = sorted(measured - described)
    if without_info:
        warnings.append(
            "These plates have measurements but no plate information, so the report "
            f"has no cell type, condition or replicate for them: {', '.join(without_info)}."
        )
    without_measurements = sorted(described - measured)
    if without_measurements:
        warnings.append(
            "These plates have plate information but no measurements, so they are "
            f"not in the report: {', '.join(without_measurements)}."
        )
    return warnings
