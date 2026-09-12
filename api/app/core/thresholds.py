"""
Shared check for wells that are running low.

The values a well reports are compared against the thresholds that are kept in
the `Threshold` model and shown on the messages page. Both the Echo import and
the `find_problems` recalculation use this check, so that a well is judged the
same way no matter which of the two looked at it.
"""


def is_below_threshold(
    current_amount: float | None,
    current_dmso: float | None,
    threshold_amount: float,
    threshold_dmso: float,
) -> bool:
    """
    Tell whether a well reported a value below one of the thresholds.
    Each value is checked on its own, so a well that reports only one of the
    two is still checked for that one. `None` means "never reported" and is
    not a problem, while zero means "empty" and is one.
    The volume is in microliter, the DMSO is a percentage.
    Example input:
    {"current_amount": 2.0, "current_dmso": 95,
     "threshold_amount": 2.5, "threshold_dmso": 80}
    Example output:
    True
    """
    if current_amount is not None and current_amount < threshold_amount:
        return True

    if current_dmso is not None and current_dmso < threshold_dmso:
        return True

    return False
