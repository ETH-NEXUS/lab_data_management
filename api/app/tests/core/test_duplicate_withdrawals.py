"""
Withdrawals that a plate mapped anew may have repeated, and the command that
lists them and removes them for a plate only when someone confirms it.

A withdrawal without target well can be a real transfer: removing it would make
LDM think a library well holds more than it does. So nothing may be removed
without --target-plate, and never more than one set per confirmation.
"""

from datetime import timedelta
from io import StringIO

from django.core.management import CommandError, call_command
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
        self.report = [(self.well_a1, 20, READING), (self.well_a2, 30, READING)]
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

    def command(self, *options):
        output = StringIO()
        call_command("remove_duplicate_withdrawals", *options, stdout=output)
        return output.getvalue()

    def without_target(self):
        return WellWithdrawal.objects.filter(target_well__isnull=True).count()

    def assertNoCandidates(self):
        self.assertIn("Plates that repeat some of them as a whole: 0", self.command())

    # Without --target-plate nothing is removed

    def test_the_list_shows_a_plate_mapped_again_and_changes_nothing(self):
        self.map_report("EXP_1", 0, self.report).delete()
        self.map_report("EXP_1", 1, self.report)

        output = self.command()

        self.assertEqual(4, WellWithdrawal.objects.count())
        self.assertIn("Plates that repeat some of them as a whole: 1", output)
        self.assertIn(
            "EXP_1 (from LIB_001): 2 withdrawals, 50.0 nL (A1, A2), repeated 1 time(s)",
            output,
        )
        self.assertIn("--target-plate <barcode> --dry-run", output)
        self.assertIn("Nothing was changed.", output)

    def test_replicate_plates_are_only_listed(self):
        # Like L3900-2_11 in 2024_snl_pruschy_radiation: four real plates with the
        # same wells, volumes and readings; the first one is deleted later
        first = self.map_report("240716MP-1_3", 0, self.report)
        self.map_report("240716MP-1_4", 0.01, self.report)
        self.map_report("240716MP-2_3", 0.02, self.report)
        self.map_report("240716MP-2_4", 0.03, self.report)
        first.delete()

        output = self.command()

        self.assertIn("Plates that repeat some of them as a whole: 3", output)
        self.assertEqual(2, self.without_target())
        self.assertEqual(8, WellWithdrawal.objects.count())

    # With --target-plate one set is removed

    def test_a_confirmed_plate_removes_one_set(self):
        self.map_report("EXP_1", 0, self.report).delete()
        self.map_report("EXP_1", 1, self.report)
        self.assertEqual(960, self.well_a1.amount)

        output = self.command("--target-plate", "EXP_1")

        self.assertEqual(0, self.without_target())
        self.assertEqual(2, WellWithdrawal.objects.count())
        self.assertEqual(980, self.well_a1.amount)
        self.assertIn("Removed 2 withdrawals that EXP_1 repeated.", output)
        self.assertNoCandidates()

    def test_a_dry_run_with_a_plate_changes_nothing(self):
        self.map_report("EXP_1", 0, self.report).delete()
        self.map_report("EXP_1", 1, self.report)

        output = self.command("--target-plate", "EXP_1", "--dry-run")

        self.assertEqual(4, WellWithdrawal.objects.count())
        self.assertIn("Would remove: EXP_1 (from LIB_001): 2 withdrawals", output)
        self.assertIn("Dry run: nothing was changed.", output)

    def test_a_report_mapped_three_times_needs_two_confirmations(self):
        self.map_report("EXP_1", 0, self.report).delete()
        self.map_report("EXP_1", 1, self.report).delete()
        self.map_report("EXP_1", 2, self.report)

        self.command("--target-plate", "EXP_1")

        self.assertEqual(2, self.without_target())
        self.assertIn("repeated 1 time(s)", self.command())

        self.command("--target-plate", "EXP_1")

        self.assertEqual(0, self.without_target())

    def test_the_newest_set_is_removed(self):
        # A real replicate from long before, and the report mapped again later
        self.map_report("EXP_0", 0, self.report).delete()
        self.map_report("EXP_1", 5, self.report).delete()
        self.map_report("EXP_1", 6, self.report)

        self.command("--target-plate", "EXP_1")

        kept = WellWithdrawal.objects.filter(target_well__isnull=True)
        self.assertEqual({self.at(0)}, {withdrawal.created_at for withdrawal in kept})

    def test_a_plate_that_repeats_only_one_well_removes_nothing(self):
        # Two different runs; by chance A1 has the same amount and reading
        self.map_report("EXP_1", 0, self.report).delete()
        self.map_report(
            "EXP_2", 1, [(self.well_a1, 20, READING), (self.well_a2, 30, OTHER_READING)]
        )

        output = self.command("--target-plate", "EXP_2")

        self.assertIn("EXP_2 does not repeat withdrawals without target well", output)
        self.assertEqual(2, self.without_target())

    def test_an_unknown_plate_is_an_error(self):
        with self.assertRaises(CommandError) as raised:
            self.command("--target-plate", "NO_SUCH_PLATE")

        self.assertIn("There is no plate NO_SUCH_PLATE.", str(raised.exception))

    # What is never a candidate

    def test_a_deleted_plate_without_new_mapping(self):
        self.map_report("EXP_1", 0, [(self.well_a1, 20, READING)]).delete()

        self.assertNoCandidates()

    def test_a_later_transfer_with_another_reading(self):
        self.map_report("EXP_1", 0, [(self.well_a1, 20, READING)]).delete()
        self.map_report("EXP_2", 1, [(self.well_a1, 20, OTHER_READING)])

        self.assertNoCandidates()

    def test_another_plate_of_the_same_run(self):
        # Echo creates all plates of a run before the withdrawals of the run
        first_plate = self.new_plate("EXP_1", 0)
        second_plate = self.new_plate("EXP_2", 0)
        self.withdraw(first_plate, 0, [(self.well_a1, 20, READING)])
        self.withdraw(second_plate, 0.001, [(self.well_a1, 20, READING)])
        first_plate.delete()

        self.assertNoCandidates()

    def test_a_second_plate_copy(self):
        # An empty well: both copies carry the reading 0 over, no PlateMapping
        empty = (0.0, 95.0)
        self.map_report(
            "LIB_001_COPY", 0, [(self.well_a1, 20, empty)], from_report=False
        ).delete()
        self.map_report(
            "LIB_001_COPY_2", 1, [(self.well_a1, 20, empty)], from_report=False
        )

        self.assertNoCandidates()

    def test_a_withdrawal_without_echo_reading(self):
        # A CSV report has no reading, so nothing tells whether it was the same transfer
        self.map_report("EXP_1", 0, [(self.well_a1, 20, None)]).delete()
        self.map_report("EXP_1", 1, [(self.well_a1, 20, None)])

        self.assertNoCandidates()

    def test_an_earlier_plate(self):
        # The withdrawal without target is the newer one, so nothing repeats it
        self.map_report("EXP_2", 0, [(self.well_a1, 20, READING)])
        self.map_report("EXP_1", 1, [(self.well_a1, 20, READING)]).delete()

        self.assertNoCandidates()
