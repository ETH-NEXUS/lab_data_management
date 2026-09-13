"""
Shared check for wells that are running low.

The values a well reports are compared against the thresholds that are kept in
the `Threshold` model and shown on the messages page. Both the Echo import and
the `find_problems` recalculation use this check, so that a well is judged the
same way no matter which of the two looked at it.
"""

VOLUME_REASON = "volume"
DMSO_REASON = "dmso"


def threshold_reasons(
    current_amount: float | None,
    current_dmso: float | None,
    threshold_amount: float,
    threshold_dmso: float,
) -> list[str]:
    """
    Name the thresholds a well is below, so that the page can say why a well
    was marked instead of only that it was.
    Each value is checked on its own, so a well that reports only one of the
    two is still checked for that one. `None` means "never reported" and is
    not a problem, while zero means "empty" and is one.
    The volume is in microliter, the DMSO is a percentage.
    Example input:
    {"current_amount": 0, "current_dmso": 0,
     "threshold_amount": 2.5, "threshold_dmso": 80}
    Example output:
    ["volume", "dmso"]
    """
    reasons = []

    if current_amount is not None and current_amount < threshold_amount:
        reasons.append(VOLUME_REASON)

    if current_dmso is not None and current_dmso < threshold_dmso:
        reasons.append(DMSO_REASON)

    return reasons


def is_below_threshold(
    current_amount: float | None,
    current_dmso: float | None,
    threshold_amount: float,
    threshold_dmso: float,
) -> bool:
    """
    Tell whether a well reported a value below one of the thresholds.
    Example input:
    {"current_amount": 2.0, "current_dmso": 95,
     "threshold_amount": 2.5, "threshold_dmso": 80}
    Example output:
    True
    """
    return bool(
        threshold_reasons(
            current_amount, current_dmso, threshold_amount, threshold_dmso
        )
    )
