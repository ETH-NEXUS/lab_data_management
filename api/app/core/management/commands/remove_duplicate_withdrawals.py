"""
Finds withdrawals that a plate mapped anew may have repeated, and removes them
for a plate only when someone confirms that it was mapped anew.

The database cannot tell a report mapped again from a real transfer to a
replicate plate with the same numbers (see core/utils/wells/duplicate_withdrawals.py),
so without --target-plate this command only lists the candidates.

Examples:
    python manage.py remove_duplicate_withdrawals
    python manage.py remove_duplicate_withdrawals --target-plate EXP_1 --dry-run
    python manage.py remove_duplicate_withdrawals --target-plate EXP_1
"""

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from core.models import Plate, PlateDetail, WellDetail, WellWithdrawal
from core.utils.wells.duplicate_withdrawals import RepeatedMapping, repeated_mappings

# How many wells of a plate the output lists before it only counts them
LISTED_WELLS = 10


class Command(BaseCommand):
    help = (
        "List plates that may repeat withdrawals of a deleted plate; remove them "
        "for one plate with --target-plate once you know it was mapped anew."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--target-plate",
            help=(
                "Barcode of a plate that was deleted and mapped again from the "
                "same report; one set of its repeated withdrawals is removed"
            ),
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="With --target-plate: only show what would be removed",
        )

    def handle(self, *args, **options):
        barcode = options["target_plate"]
        if barcode is None:
            self.list_candidates()
            return

        if not Plate.objects.filter(barcode=barcode).exists():
            raise CommandError(f"There is no plate {barcode}.")
        repeated = repeated_mappings(target_barcode=barcode)
        if not repeated:
            self.stdout.write(
                f"{barcode} does not repeat withdrawals without target well. "
                "Nothing was changed."
            )
            return

        label = "Would remove" if options["dry_run"] else "Removing"
        to_remove = []
        for mapping in repeated:
            self.stdout.write(f"{label}: {self.describe(mapping)}")
            to_remove.extend(mapping.one_set)
        if options["dry_run"]:
            self.stdout.write("Dry run: nothing was changed.")
            return

        with transaction.atomic():
            WellWithdrawal.objects.filter(
                id__in=[withdrawal.id for withdrawal in to_remove]
            ).delete()
        # The remaining amount of a well in the views subtracts all withdrawals
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        self.stdout.write(
            self.style.SUCCESS(
                f"Removed {len(to_remove)} withdrawals that {barcode} repeated."
            )
        )

    def list_candidates(self) -> None:
        """Prints the plates that may repeat withdrawals; changes nothing."""
        without_target = WellWithdrawal.objects.filter(target_well__isnull=True)
        repeated = repeated_mappings()
        self.stdout.write(f"Withdrawals without target well: {without_target.count()}")
        self.stdout.write(
            f"Plates that repeat some of them as a whole: {len(repeated)}"
        )
        for mapping in repeated:
            self.stdout.write(f"  {self.describe(mapping)}")
        if repeated:
            self.stdout.write(
                "Such a plate is either the same report mapped again (the "
                "withdrawals are counted twice) or a real transfer with the same "
                "numbers, e.g. a replicate plate. The database cannot tell which. "
                "Only if you know that a plate was deleted and mapped again, run:\n"
                "  python manage.py remove_duplicate_withdrawals "
                "--target-plate <barcode> --dry-run"
            )
        self.stdout.write("Nothing was changed.")

    def describe(self, mapping: RepeatedMapping) -> str:
        """
        One line per plate, e.g. "EXP_1 (from LIB_001): 2 withdrawals, 50.0 nL
        (A1, A2), repeated 1 time(s)".
        """
        total = sum(withdrawal.amount for withdrawal in mapping.one_set)
        wells = [withdrawal.well.hr_position for withdrawal in mapping.one_set]
        listed = ", ".join(wells[:LISTED_WELLS])
        if len(wells) > LISTED_WELLS:
            listed += f" and {len(wells) - LISTED_WELLS} more"
        return (
            f"{mapping.target_plate.barcode} (from {mapping.source_plate.barcode}): "
            f"{len(mapping.one_set)} withdrawals, {total} nL ({listed}), "
            f"repeated {mapping.full_sets} time(s)"
        )
