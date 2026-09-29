"""
The control wells of a report. The R reports know two well types as controls:
"P" (positive) and "N" (negative). In LDM the controls may be named differently,
e.g. "P1" or "N2", so the user chooses which well types are the controls, and
they are renamed to "P" and "N" in the input files.
"""

import pandas as pd
from django.core.management.base import CommandError


def rename_controls(frame: pd.DataFrame, positive: str, negative: str) -> pd.DataFrame:
    """
    The rows with the "control" column (the well type) as the report reads it:
    the chosen controls become "P" and "N". Another well type that is called "P"
    or "N" but was not chosen gets " (not chosen)", so R does not take it as a control.

    positive = "P1", negative = "N1":
    "P1" -> "P", "N1" -> "N", "C" -> "C", "P" -> "P (not chosen)"
    """
    renamed = frame.copy()
    new_names = {positive: "P", negative: "N"}

    def report_name(well_type: str) -> str:
        if well_type in new_names:
            return new_names[well_type]
        if well_type in ["P", "N"]:
            return f"{well_type} (not chosen)"
        return well_type

    renamed["control"] = renamed["control"].map(report_name)
    return renamed


def check_chosen_controls(
    main_info: pd.DataFrame, positive: str, negative: str
) -> None:
    """The chosen controls must be two different well types of this measurement."""
    well_types = sorted(set(main_info["control"]))
    # e.g. '"C", "N1", "P1"', in quotes so a space at the end is seen
    well_types_text = ", ".join('"' + well_type + '"' for well_type in well_types)
    if positive == negative:
        raise CommandError(
            f'The positive and the negative control are both "{positive}". Choose two '
            "different well types in the analysis window."
        )
    for name, chosen in [("positive", positive), ("negative", negative)]:
        if chosen not in well_types:
            raise CommandError(
                f'The {name} control "{chosen}" is not a well type of this measurement. '
                f"Its well types are: {well_types_text}. "
                "Choose the control wells in the analysis window."
            )


def control_warnings(
    main_info: pd.DataFrame, positive: str, negative: str
) -> list[str]:
    """
    The report normalizes every plate by its positive and negative control wells;
    a plate without them is left out of the results (seen on experiment 105:
    plate 20250513SP_29 has no P wells and no results). `main_info` has the
    renamed controls already ("P" and "N").

    Returned data example:
    ['Plate SP_2 has no negative control wells ("N1"), so the report ...']
    """
    names = {
        "P": f'positive control wells ("{positive}")',
        "N": f'negative control wells ("{negative}")',
    }
    warnings = []
    plates = main_info.groupby("plate")["control"]
    for plate, controls in plates:
        missing = [
            names[control] for control in ["P", "N"] if control not in set(controls)
        ]
        if missing:
            warnings.append(
                f"Plate {plate} has no {' and no '.join(missing)}, so the report "
                "cannot normalize it and leaves it out of the results."
            )
    # Without a single plate to normalize, R stops with an error that does not say why
    if len(warnings) == len(plates):
        raise CommandError(
            f'No plate of this measurement has both the positive ("{positive}") and '
            f'the negative ("{negative}") control wells, so the report cannot '
            "normalize any plate. Check the controls chosen in the analysis window "
            "and the well types of the plate layout.\n" + "\n".join(warnings)
        )
    return warnings
