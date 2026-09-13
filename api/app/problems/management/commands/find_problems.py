from django.core.management import BaseCommand
from django.db import transaction
from django.db.models import Count, Prefetch

from core.models import Plate, Threshold, Well, WellWithdrawal
from core.thresholds import is_below_threshold


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("what", help="Which problem to solve")  # mark_empty_wells

    def mark_empty_wells(self):
        threshold = Threshold.current()

        # The newest withdrawal carries the current state of the well.
        # We order the prefetch so that the first entry is the newest one,
        # the same way Well.current_info picks its withdrawal.
        newest_withdrawals_first = WellWithdrawal.objects.order_by("-created_at")
        library_wells = (
            Well.objects.filter(plate__library__isnull=False)
            .annotate(withdrawals_count=Count("withdrawals"))
            .filter(withdrawals_count__gt=0)
            .select_related("plate__dimension")
            .prefetch_related(
                Prefetch("withdrawals", queryset=newest_withdrawals_first)
            )
        )
        print(f"Number of library wells: {library_wells.count()}")

        wells_to_update = []
        problematic_plate_ids = set()

        for well in library_wells:
            withdrawals = list(well.withdrawals.all())
            last_withdrawal = withdrawals[0] if withdrawals else None

            if last_withdrawal is not None and is_below_threshold(
                last_withdrawal.current_amount,
                last_withdrawal.current_dmso,
                threshold.amount,
                threshold.dmso,
            ):
                problematic_plate_ids.add(well.plate_id)
                print(f"Marking well {well.hr_position} as empty")
                print(f"Amount: {last_withdrawal.current_amount}")
                print(f"DMSO: {last_withdrawal.current_dmso}")
                if well.status != "empty":
                    well.status = "empty"
                    wells_to_update.append(well)
            else:
                if well.status == "empty":
                    well.status = None
                    wells_to_update.append(well)

        # The status of a plate depends on all of its wells, so it can only be
        # decided after every well has been looked at.
        plates_to_flag = Plate.objects.filter(id__in=problematic_plate_ids).exclude(
            status="empty_wells"
        )
        plates_to_unflag = Plate.objects.filter(
            library__isnull=False, status="empty_wells"
        ).exclude(id__in=problematic_plate_ids)

        print(f"Number of wells to update: {len(wells_to_update)}")
        print(f"Number of plates to flag: {plates_to_flag.count()}")
        print(f"Number of plates to unflag: {plates_to_unflag.count()}")

        with transaction.atomic():
            Well.objects.bulk_update(wells_to_update, ["status"], batch_size=1000)
            plates_to_flag.update(status="empty_wells")
            plates_to_unflag.update(status=None)

    def handle(self, *args, **options):
        if options["what"] == "mark_empty_wells":
            self.mark_empty_wells()
