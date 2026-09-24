"""
Deleting a plate together with the withdrawals that were transferred into it.

The lab deletes a plate only to correct a mapping that went wrong in LDM, and then
maps the same Echo report again. Without this, the withdrawals from the library wells
would stay (their `target_well` becomes empty) and the new mapping would add them a
second time, so the library wells would lose the volume twice.

Deleting an experiment or a project does not go through here: its plates keep the
old behaviour, and the withdrawals into them stay as a record of real transfers.
"""

from django.db import transaction

from core.models import Plate, WellWithdrawal


def delete_plate_with_withdrawals(plate: Plate) -> None:
    """Deletes the plate and the withdrawals into its wells, all or nothing."""
    with transaction.atomic():
        WellWithdrawal.objects.filter(target_well__plate=plate).delete()
        plate.delete()
