"""
Withdrawals counted twice after an experiment plate was deleted and the same
Echo report was mapped again, and the command that removes them.
"""

from datetime import timedelta
from io import StringIO

from django.core.management import call_command
from django.test import TestCase
from django.utils import timezone

from compoundlib.models import Compound
from core.models import Plate, PlateDimension, Well, WellCompound, WellWithdrawal


class DuplicateWithdrawalsTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        library_plate = Plate.objects.create(
            barcode="LIB_001", dimension=self.dimension
        )
        self.library_well = Well.objects.create(plate=library_plate, position=0)
        WellCompound.objects.create(
            well=self.library_well,
            compound=Compound.objects.create(name="Compound A"),
            amount=1000,
        )
        self.start = timezone.now() - timedelta(days=10)

    def transfer(self, barcode, days_later, amount=20, reading=(900.0, 95.0)):
        """A transfer from the library well to a new experiment plate."""
        plate = Plate.objects.create(barcode=barcode, dimension=self.dimension)
        target_well = Well.objects.create(plate=plate, position=0)
        withdrawal = WellWithdrawal.objects.create(
            well=self.library_well,
            target_well=target_well,
            amount=amount,
            current_amount=reading[0] if reading else None,
            current_dmso=reading[1] if reading else None,
        )
        # The order of the withdrawals is what counts, so it is set explicitly
        WellWithdrawal.objects.filter(id=withdrawal.id).update(
            created_at=self.start + timedelta(days=days_later)
        )
        return plate

    def remove(self, *options):
        output = StringIO()
        call_command("remove_duplicate_withdrawals", *options, stdout=output)
        return output.getvalue()

    def test_the_same_report_mapped_again_is_removed(self):
        self.transfer("EXP_1", days_later=0).delete()
        self.transfer("EXP_1", days_later=1)
        self.assertEqual(960, self.library_well.amount)

        output = self.remove()

        self.assertEqual(1, WellWithdrawal.objects.count())
        self.assertIsNotNone(WellWithdrawal.objects.get().target_well)
        self.assertEqual(980, self.library_well.amount)
        self.assertIn("LIB_001: 1 duplicates, 20.0 nL (A1)", output)
        self.assertIn("Removed 1 duplicate withdrawals.", output)

    def test_a_dry_run_changes_nothing(self):
        self.transfer("EXP_1", days_later=0).delete()
        self.transfer("EXP_1", days_later=1)

        output = self.remove("--dry-run")

        self.assertEqual(2, WellWithdrawal.objects.count())
        self.assertIn("Of these, duplicates: 1", output)
        self.assertIn("Nothing was changed.", output)

    def test_a_later_transfer_with_another_reading_is_a_real_one(self):
        # Same standard volume, but the well held less the second time
        self.transfer("EXP_1", days_later=0, reading=(900.0, 95.0)).delete()
        self.transfer("EXP_2", days_later=1, reading=(880.0, 95.0))

        self.remove()

        self.assertEqual(2, WellWithdrawal.objects.count())

    def test_a_deleted_plate_without_new_mapping_keeps_its_withdrawal(self):
        # The liquid was transferred, the plate was deleted later for another reason
        self.transfer("EXP_1", days_later=0).delete()

        output = self.remove()

        self.assertEqual(1, WellWithdrawal.objects.count())
        self.assertIn("Withdrawals without target well: 1", output)
        self.assertIn("Of these, duplicates: 0", output)

    def test_a_withdrawal_without_echo_reading_is_kept(self):
        # A CSV report has no reading, so nothing tells whether it was the same transfer
        self.transfer("EXP_1", days_later=0, reading=None).delete()
        self.transfer("EXP_1", days_later=1, reading=None)

        self.remove()

        self.assertEqual(2, WellWithdrawal.objects.count())

    def test_an_earlier_withdrawal_with_target_is_not_a_repetition(self):
        # The withdrawal without target is the newer one, so it is not repeated
        self.transfer("EXP_2", days_later=0)
        self.transfer("EXP_1", days_later=1).delete()

        self.remove()

        self.assertEqual(2, WellWithdrawal.objects.count())
