"""
Volume units used in this project.

Two different units live next to each other, because the Echo transfer files
report them that way:

- transferred and stored volumes (`WellCompound.amount`, `WellWithdrawal.amount`,
  the plate volumes in the admin) are in nanoliter (nL),
- the current fill level of a well (`WellWithdrawal.current_amount`) is in
  microliter (uL). This is the unit of the "Current Fluid Volume" column of the
  Echo files and the unit of the volume threshold on the messages page.

Whenever a stored volume is used as a fill level, it has to be converted here.
"""

NANOLITERS_PER_MICROLITER = 1000


def nanoliter_to_microliter(volume_in_nanoliter: float) -> float:
    """
    Convert a transferred volume into the unit used for fill levels.
    Example input:
    {"volume_in_nanoliter": 2500}
    Example output:
    2.5
    """
    return volume_in_nanoliter / NANOLITERS_PER_MICROLITER
