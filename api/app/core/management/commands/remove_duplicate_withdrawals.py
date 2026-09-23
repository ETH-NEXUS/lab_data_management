"""
Removes the withdrawals that a new mapping of the same Echo report repeated.

See core/utils/wells/duplicate_withdrawals.py for which withdrawals count as
duplicates: only those that a plate mapped anew repeats as a whole. Every other
withdrawal without target well stays: it may be a real transfer to a plate that
was deleted later.

Example:
    python manage.py remove_duplicate_withdrawals --dry-run
"""

from django.core.management.base import BaseCommand
from django.db import transaction

from core.models import PlateDetail, WellDetail, WellWithdrawal
from core.utils.wells.duplicate_withdrawals import (
    RepeatedMapping,
    repeated_mappings,
)

# How many wells of a plate the report lists before it only counts them
LISTED_WELLS = 10


class Command(BaseCommand):
    help = "Remove withdrawals that a new mapping of the same Echo report repeated."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Only report what would be removed, change nothing",
        )

    def handle(self, *args, **options):
        repeated = repeated_mappings()
        duplicates = []
        for mapping in repeated:
            duplicates.extend(mapping.duplicates)
        without_target = WellWithdrawal.objects.filter(target_well__isnull=True)
        self.report(repeated, len(duplicates), without_target.count())

        if options["dry_run"] or not duplicates:
            self.stdout.write("Nothing was changed.")
            return

        with transaction.atomic():
            WellWithdrawal.objects.filter(
                id__in=[withdrawal.id for withdrawal in duplicates]
            ).delete()
        # The remaining amount of a well in the views subtracts all withdrawals
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        self.stdout.write(
            self.style.SUCCESS(f"Removed {len(duplicates)} duplicate withdrawals.")
        )

    def report(
        self,
        repeated: list[RepeatedMapping],
        duplicate_count: int,
        without_target_count: int,
    ) -> None:
        """
        Prints the duplicates by repeated mapping, e.g.
        "LIB_001 -> EXP_1 (mapped anew): 3 duplicates, 60.0 nL (A1, A2, B7)".
        """
        self.stdout.write(f"Withdrawals without target well: {without_target_count}")
        self.stdout.write(f"Of these, duplicates: {duplicate_count}")

        for mapping in repeated:
            total = sum(withdrawal.amount for withdrawal in mapping.duplicates)
            wells = [withdrawal.well.hr_position for withdrawal in mapping.duplicates]
            listed = ", ".join(wells[:LISTED_WELLS])
            if len(wells) > LISTED_WELLS:
                listed += f" and {len(wells) - LISTED_WELLS} more"
            self.stdout.write(
                f"{mapping.source_plate.barcode} -> {mapping.target_plate.barcode} "
                f"(mapped anew): {len(mapping.duplicates)} duplicates, {total} nL "
                f"({listed})"
            )
