"""
Tests for the `find_problems mark_empty_wells` recalculation.

The command decides two things: which wells are running low, and which plates
therefore have to show up on the messages page. The plate part is what used to
go wrong - a single well that was fine unflagged the whole plate again.
"""

from datetime import timedelta

from django.core.management import CommandError, call_command
from django.test import SimpleTestCase, TestCase
from django.utils import timezone

from compoundlib.models import Compound, CompoundLibrary
from core.models import (
    Plate,
    PlateDimension,
    Threshold,
    Well,
    WellCompound,
    WellWithdrawal,
)


class MarkEmptyWellsTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        Threshold.objects.create(amount=2.5, dmso=80)
        self.dimension = PlateDimension.objects.create(name="dim_4x2", cols=4, rows=2)
        self.library = CompoundLibrary.objects.create(name="Test library")
        self.plate = Plate.objects.create(
            barcode="LIB_PLATE", dimension=self.dimension, library=self.library
        )
        self.compound = Compound.objects.create(name="Test compound")

    def add_well(self, position, reports, status=None):
        """
        Adds a well and the values it reported, oldest first.
        Reports example: [(9.5, 95), (2.0, 95)]
        """
        well = Well.objects.create(plate=self.plate, position=position, status=status)
        WellCompound.objects.create(well=well, compound=self.compound, amount=1000)
        now = timezone.now()
        for index, (current_amount, current_dmso) in enumerate(reports):
            WellWithdrawal.objects.create(
                well=well,
                target_well=None,
                amount=30,
                current_amount=current_amount,
                current_dmso=current_dmso,
                created_at=now - timedelta(days=len(reports) - index),
            )
        return well

    def recalculate(self):
        call_command("find_problems", "mark_empty_wells", verbosity=0)
        self.plate.refresh_from_db()

    def well_status(self, well):
        well.refresh_from_db()
        return well.status

    def test_a_low_volume_is_marked(self):
        well = self.add_well(0, [(2.0, 95)])
        self.recalculate()
        self.assertEqual("empty", self.well_status(well))
        self.assertEqual("empty_wells", self.plate.status)

    def test_a_volume_of_zero_is_marked(self):
        well = self.add_well(1, [(0, 95)])
        self.recalculate()
        self.assertEqual("empty", self.well_status(well))

    def test_a_low_dmso_is_marked(self):
        well = self.add_well(2, [(9.5, 70)])
        self.recalculate()
        self.assertEqual("empty", self.well_status(well))

    def test_a_well_that_is_fine_is_not_marked(self):
        well = self.add_well(3, [(9.5, 95)])
        self.recalculate()
        self.assertIsNone(self.well_status(well))
        self.assertIsNone(self.plate.status)

    def test_the_newest_report_decides(self):
        """The well was refilled after the report that was below the threshold."""
        well = self.add_well(4, [(2.0, 95), (9.5, 95)])
        self.recalculate()
        self.assertIsNone(self.well_status(well))

    def test_a_well_that_is_fine_again_loses_its_mark(self):
        well = self.add_well(5, [(9.5, 95)], status="empty")
        self.recalculate()
        self.assertIsNone(self.well_status(well))

    def test_one_low_well_keeps_the_plate_flagged_next_to_healthy_ones(self):
        low_well = self.add_well(0, [(2.0, 95)])
        for position in (1, 2, 3):
            self.add_well(position, [(9.5, 95)])

        self.recalculate()
        self.assertEqual("empty_wells", self.plate.status)

        # Running it again must not unflag the plate: that is what used to hide
        # plates with only a few problematic wells from the messages page.
        self.recalculate()
        self.assertEqual("empty_wells", self.plate.status)
        self.assertEqual("empty", self.well_status(low_well))

    def test_a_plate_without_problems_is_unflagged(self):
        Plate.objects.filter(id=self.plate.id).update(status="empty_wells")
        self.add_well(0, [(9.5, 95)])
        self.recalculate()
        self.assertIsNone(self.plate.status)

    def test_wells_of_a_plate_without_a_library_are_left_alone(self):
        plate = Plate.objects.create(barcode="NO_LIBRARY", dimension=self.dimension)
        well = Well.objects.create(plate=plate, position=0)
        WellWithdrawal.objects.create(
            well=well, target_well=None, amount=30, current_amount=0, current_dmso=0
        )
        self.recalculate()
        self.assertIsNone(self.well_status(well))

    def test_a_status_set_by_hand_is_kept(self):
        Plate.objects.filter(id=self.plate.id).update(status="disposed")
        low_well = self.add_well(0, [(0, 0)])
        self.recalculate()
        self.assertEqual("disposed", self.plate.status)
        self.assertEqual("empty", self.well_status(low_well))

    def test_a_plate_with_an_empty_status_is_flagged(self):
        Plate.objects.filter(id=self.plate.id).update(status="")
        self.add_well(0, [(2.0, 95)])
        self.recalculate()
        self.assertEqual("empty_wells", self.plate.status)


class UnknownProblemTest(SimpleTestCase):
    def test_an_unknown_problem_is_an_error(self):
        # A typo used to do nothing and report nothing
        with self.assertRaises(CommandError) as raised:
            call_command("find_problems", "mark_empty_well")

        self.assertIn("invalid choice: 'mark_empty_well'", str(raised.exception))
