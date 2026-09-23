"""
Withdrawals that were counted twice.

Deleting the plates of an experiment and mapping the same Echo report again is a
common way to repair a mapping. The withdrawals of the library wells are kept
when their target plate is deleted (`target_well` becomes empty), so the second
mapping adds the same withdrawals again.

A withdrawal without target well can also be a real one: the liquid was
transferred, and the experiment plate was deleted later for another reason.
So only a clear repetition counts as duplicate (see `duplicate_withdrawals`).
"""

from django.db.models import Exists, OuterRef, QuerySet

from core.models import WellWithdrawal


def duplicate_withdrawals() -> QuerySet:
    """
    The withdrawals without target well that a later withdrawal of the same
    well repeats: same amount and same Echo reading (current_amount and
    current_dmso of the source well at the transfer).

    The reading is what tells two transfers apart: the amounts are standard
    volumes (10, 20, 30 nL ...), but the volume left in the well is lower after
    every real transfer. A withdrawal without reading is never a duplicate here,
    because nothing tells whether it was the same transfer.
    """
    later_repetition = WellWithdrawal.objects.filter(
        well=OuterRef("well"),
        target_well__isnull=False,
        created_at__gt=OuterRef("created_at"),
        amount=OuterRef("amount"),
        current_amount=OuterRef("current_amount"),
        current_dmso=OuterRef("current_dmso"),
    )
    return (
        WellWithdrawal.objects.filter(target_well__isnull=True)
        .filter(Exists(later_repetition))
        .select_related("well__plate__dimension")
    )
