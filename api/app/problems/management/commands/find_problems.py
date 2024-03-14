from django.core.management import BaseCommand
from django.db.models import F, Count

from core.models import Threshold, Well, Plate


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("what", help="Which problem to solve")  # mark_empty_wells

    def mark_empty_wells(self):
        threshold = Threshold.objects.first()
        if not threshold:
            return

        library_wells = (
            Well.objects.filter(plate__library__isnull=False)
            .annotate(withdrawals_count=Count("withdrawals"))
            .filter(withdrawals_count__gt=0)
            .select_related("plate")
            .prefetch_related("withdrawals")
        )
        print(f"Number of library wells: {library_wells.count()}")

        wells_to_update = []
        plates_to_update = set()

        for well in library_wells:
            last_withdrawal = well.withdrawals.last()
            if (
                not last_withdrawal
                or not last_withdrawal.current_amount
                or not last_withdrawal.current_dmso
            ):
                continue

            plate = well.plate
            well_changed = False
            plate_changed = False

            if (
                last_withdrawal.current_amount < threshold.amount
                or last_withdrawal.current_dmso < threshold.dmso
            ):
                print(f"Marking well {well.hr_position} as empty")
                print(f"Amount: {last_withdrawal.current_amount}")
                print(f"DMSO: {last_withdrawal.current_dmso}")
                if well.status != "empty":
                    well.status = "empty"
                    well_changed = True
                if plate.status != "empty_wells":
                    plate.status = "empty_wells"
                    plate_changed = True
            else:
                if well.status == "empty":
                    well.status = None
                    well_changed = True
                if plate.status == "empty_wells":
                    plate.status = None
                    plate_changed = True

            if well_changed:
                wells_to_update.append(well)
            if plate_changed:
                plates_to_update.add(plate)
        print(f"Number of wells to update: {len(wells_to_update)}")
        print(f"Number of plates to update: {len(plates_to_update)}")
        Well.objects.bulk_update(wells_to_update, ["status"])
        Plate.objects.bulk_update(plates_to_update, ["status"])

    def handle(self, *args, **options):
        if options["what"] == "mark_empty_wells":
            self.mark_empty_wells()
