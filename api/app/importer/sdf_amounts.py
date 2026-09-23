"""
The amounts of an SDF library in nanoliter, the unit LDM stores well amounts
in (the Echo withdrawals are in nanoliter as well).

Only the "Vol_Copy…" columns are known to be volumes: the lab confirmed in
September 2026 that they are in microliter (e.g. "6" or "24.0"). Other amount
columns are not used, because they can mean something else: PLATE_AMOUNT1 of
a vendor file is close to the molecular weight, so probably a mass in µg.
"""

import math
from collections import Counter

VOLUME_COLUMN_PREFIX = "Vol_Copy"
NANOLITER_PER_MICROLITER = 1000


def is_volume_column(column_name: str) -> bool:
    """True for the columns with a volume in µL, e.g. "Vol_Copy1"."""
    return column_name.startswith(VOLUME_COLUMN_PREFIX)


def amount_in_nanoliter(value: str | float | None) -> float | None:
    """
    A volume in µL from the SDF file, in nL:
    "6" -> 6000.0, "24.0" -> 24000.0.
    An empty value is a copy without this compound -> 0.0.
    A value that is not a number, e.g. "<24" (less than 24 µL) -> None.
    """
    # RDKit reads the properties as text; a property that a record does not
    # have is NaN.
    if value is None:
        return 0.0
    if isinstance(value, float) and math.isnan(value):
        return 0.0
    text = str(value).strip()
    if text == "":
        return 0.0
    try:
        microliter = float(text)
    except ValueError:
        return None
    return microliter * NANOLITER_PER_MICROLITER


def unknown_unit_warning(column_name: str) -> str:
    return (
        f"The amounts in column {column_name} are not stored (set to 0), because "
        f"only the {VOLUME_COLUMN_PREFIX}… columns are known to be volumes in µL."
    )


def not_a_number_warning(column_name: str, values: Counter) -> str:
    """
    One warning for all wells of a column without an exact volume,
    e.g. values = Counter({"<24": 2486}).
    """
    described_values = ", ".join(
        f"'{value}' in {count} wells" for value, count in values.items()
    )
    return (
        f"Column {column_name} has no exact volume ({described_values}), "
        f"so these amounts are set to 0."
    )
