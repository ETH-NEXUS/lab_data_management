"""
Withdrawals counted twice after an experiment plate was deleted and the same
Echo report was mapped again, and the command that removes them.

Every case where a withdrawal without target well is a real transfer must keep
it: removing a real withdrawal would make LDM think a library well holds more
than it does.
"""

from datetime import timedelta
from io import StringIO

from django.core.management import call_command
from django.test import TestCase
from django.utils import timezone

from compoundlib.models import Compound
from core.models import (
    Plate,
    PlateDimension,
    PlateMapping,
    Well,
    WellCompound,
    WellWithdrawal,
)

READING = (900.0, 95.0)
OTHER_READING = (880.0, 95.0)


class DuplicateWithdrawalsTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        self.library_plate = Plate.objects.create(
            barcode="LIB_001", dimension=self.dimension
        )
        self.well_a1 = self.library_well(0)
        self.well_a2 = self.library_well(1)
        self.start = timezone.now() - timedelta(days=10)

    def library_well(self, position):
        well = Well.objects.create(plate=self.library_plate, position=position)
        WellCompound.objects.create(
            well=well,
            compound=Compound.objects.create(name=f"Compound {position}"),
            amount=1000,
        )
        return well

    def at(self, days_later):
        return self.start + timedelta(days=days_later)

    def new_plate(self, barcode, days_later, from_report=True):
        """
        A target plate like a mapping creates it: a minute before its
        withdrawals, with a PlateMapping if it came from a report (a plate
        copy has none).
        """
        plate = Plate.objects.create(barcode=barcode, dimension=self.dimension)
        Plate.objects.filter(id=plate.id).update(
            created_at=self.at(days_later) - timedelta(minutes=1)
        )
        if from_report:
            PlateMapping.objects.create(
                source_plate=self.library_plate, target_plate=plate
            )
        plate.refresh_from_db()
        return plate

    def withdraw(self, plate, days_later, transfers):
        """
        Transfers into `plate`, one target well each, e.g.
        [(self.well_a1, 20, READING), (self.well_a2, 30, READING)].
        """
        for position, (source_well, amount, reading) in enumerate(transfers):
            target_well = Well.objects.create(plate=plate, position=position)
            withdrawal = WellWithdrawal.objects.create(
                well=source_well,
                target_well=target_well,
                amount=amount,
                current_amount=reading[0] if reading else None,
                current_dmso=reading[1] if reading else None,
            )
            # The order in time is what counts, so it is set explicitly
            WellWithdrawal.objects.filter(id=withdrawal.id).update(
                created_at=self.at(days_later)
            )
        return plate

    def map_report(self, barcode, days_later, transfers, from_report=True):
        plate = self.new_plate(barcode, days_later, from_report)
        return self.withdraw(plate, days_later, transfers)

    def remove(self, *options):
        output = StringIO()
        call_command("remove_duplicate_withdrawals", *options, stdout=output)
        return output.getvalue()

    def without_target(self):
        return WellWithdrawal.objects.filter(target_well__isnull=True).count()

    # A report mapped again

    def test_the_same_report_mapped_again_is_removed(self):
        report = [(self.well_a1, 20, READING), (self.well_a2, 30, READING)]
        self.map_report("EXP_1", 0, report).delete()
        self.map_report("EXP_1", 1, report)
        self.assertEqual(960, self.well_a1.amount)

        output = self.remove()

        self.assertEqual(2, WellWithdrawal.objects.count())
        self.assertEqual(0, self.without_target())
        self.assertEqual(980, self.well_a1.amount)
        self.assertIn(
            "LIB_001 -> EXP_1 (mapped anew): 2 duplicates, 50.0 nL (A1, A2)", output
        )
        self.assertIn("Removed 2 duplicate withdrawals.", output)

    def test_a_report_mapped_three_times_keeps_one_set(self):
        report = [(self.well_a1, 20, READING), (self.well_a2, 30, READING)]
        self.map_report("EXP_1", 0, report).delete()
        self.map_report("EXP_1", 1, report).delete()
        self.map_report("EXP_1", 2, report)

        self.remove()

        self.assertEqual(2, WellWithdrawal.objects.count())
        self.assertEqual(0, self.without_target())

    def test_a_second_run_finds_nothing_more(self):
        report = [(self.well_a1, 20, READING)]
        self.map_report("EXP_1", 0, report).delete()
        self.map_report("EXP_1", 1, report)
        self.remove()

        output = self.remove()

        self.assertIn("Of these, duplicates: 0", output)
        self.assertEqual(1, WellWithdrawal.objects.count())

    def test_a_dry_run_changes_nothing(self):
        report = [(self.well_a1, 20, READING)]
        self.map_report("EXP_1", 0, report).delete()
        self.map_report("EXP_1", 1, report)

        output = self.remove("--dry-run")

        self.assertEqual(2, WellWithdrawal.objects.count())
        self.assertIn("Of these, duplicates: 1", output)
        self.assertIn("Nothing was changed.", output)

    # Real transfers that must stay

    def test_a_deleted_plate_without_new_mapping_keeps_its_withdrawals(self):
        self.map_report("EXP_1", 0, [(self.well_a1, 20, READING)]).delete()

        output = self.remove()

        self.assertEqual(1, self.without_target())
        self.assertIn("Withdrawals without target well: 1", output)
        self.assertIn("Of these, duplicates: 0", output)

    def test_a_later_run_that_agrees_in_one_well_only_is_a_real_one(self):
        # Two different runs; by chance A1 has the same amount and reading
        self.map_report(
            "EXP_1", 0, [(self.well_a1, 20, READING), (self.well_a2, 30, READING)]
        ).delete()
        self.map_report(
            "EXP_2", 1, [(self.well_a1, 20, READING), (self.well_a2, 30, OTHER_READING)]
        )

        self.remove()

        self.assertEqual(2, self.without_target())

    def test_a_later_transfer_with_another_reading_is_a_real_one(self):
        self.map_report("EXP_1", 0, [(self.well_a1, 20, READING)]).delete()
        self.map_report("EXP_2", 1, [(self.well_a1, 20, OTHER_READING)])

        self.remove()

        self.assertEqual(1, self.without_target())

    def test_another_plate_of_the_same_echo_run_is_not_a_repetition(self):
        # One run stamps the same volume into two plates with one reading; Echo
        # creates both plates before the withdrawals of the run
        first_plate = self.new_plate("EXP_1", 0)
        second_plate = self.new_plate("EXP_2", 0)
        self.withdraw(first_plate, 0, [(self.well_a1, 20, READING)])
        self.withdraw(second_plate, 0.001, [(self.well_a1, 20, READING)])
        # The first plate is deleted later, the transfer was real
        first_plate.delete()

        self.remove()

        self.assertEqual(1, self.without_target())

    def test_a_second_plate_copy_is_not_a_repetition(self):
        # An empty well: both copies carry the reading 0 over, the copies were real
        empty = (0.0, 95.0)
        self.map_report(
            "LIB_001_COPY", 0, [(self.well_a1, 20, empty)], from_report=False
        ).delete()
        self.map_report(
            "LIB_001_COPY_2", 1, [(self.well_a1, 20, empty)], from_report=False
        )

        self.remove()

        self.assertEqual(1, self.without_target())

    def test_a_withdrawal_without_echo_reading_is_kept(self):
        # A CSV report has no reading, so nothing tells whether it was the same transfer
        self.map_report("EXP_1", 0, [(self.well_a1, 20, None)]).delete()
        self.map_report("EXP_1", 1, [(self.well_a1, 20, None)])

        self.remove()

        self.assertEqual(1, self.without_target())

    def test_an_earlier_plate_is_not_a_repetition(self):
        # The withdrawal without target is the newer one, so nothing repeats it
        self.map_report("EXP_2", 0, [(self.well_a1, 20, READING)])
        self.map_report("EXP_1", 1, [(self.well_a1, 20, READING)]).delete()

        self.remove()

        self.assertEqual(1, self.without_target())

    def test_a_real_transfer_that_agrees_with_a_repeated_one_stays(self):
        # EXP_1 was mapped again as EXP_1 (a duplicate); EXP_9 was a real run that
        # was deleted later and agrees with it in A1 by chance
        report = [(self.well_a1, 20, READING), (self.well_a2, 30, READING)]
        self.map_report("EXP_1", 0, report).delete()
        self.map_report("EXP_9", 0.5, [(self.well_a1, 20, READING)]).delete()
        self.map_report("EXP_1", 1, report)

        self.remove()

        # One set is gone; one withdrawal of A1 without target is left, and with
        # the same amount it keeps the real volume of EXP_9
        self.assertEqual(1, self.without_target())
        self.assertEqual(3, WellWithdrawal.objects.count())
        self.assertEqual(960, self.well_a1.amount)
