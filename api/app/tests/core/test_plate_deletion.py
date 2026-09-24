"""
Deleting a plate also deletes the withdrawals into it, so that mapping the same
Echo report again does not take the volume from the library wells twice.
Deleting an experiment keeps them, as before.
"""

from django.contrib.admin.sites import AdminSite
from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.test import APITestCase

from compoundlib.models import Compound, CompoundLibrary
from core.admin import PlateAdmin
from core.models import (
    Experiment,
    Plate,
    PlateDimension,
    Project,
    Well,
    WellCompound,
    WellWithdrawal,
)


class PlateDeletionTest(APITestCase):
    fixtures = ("well_types",)

    def setUp(self):
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        library = CompoundLibrary.objects.create(name="Library")
        library_plate = Plate.objects.create(
            barcode="LIB_001", dimension=self.dimension, library=library
        )
        self.library_well = Well.objects.create(plate=library_plate, position=0)
        WellCompound.objects.create(
            well=self.library_well,
            compound=Compound.objects.create(name="Aspirin"),
            amount=1000,
        )
        project = Project.objects.create(name="Project")
        self.experiment = Experiment.objects.create(name="Experiment", project=project)
        self.user = get_user_model().objects.create(username="tester")
        self.client.force_authenticate(user=self.user)

    def mapped_plate(self, barcode, amount=20):
        """A plate of the experiment with one transfer from the library well."""
        plate = Plate.objects.create(
            barcode=barcode, dimension=self.dimension, experiment=self.experiment
        )
        target_well = Well.objects.create(plate=plate, position=0)
        WellWithdrawal.objects.create(
            well=self.library_well, target_well=target_well, amount=amount
        )
        return plate

    def test_deleting_a_plate_through_the_api_deletes_the_withdrawals_into_it(self):
        plate = self.mapped_plate("EXP_1")
        self.mapped_plate("EXP_2", amount=30)

        response = self.client.delete(f"/api/plates/{plate.id}/")

        self.assertEqual(status.HTTP_204_NO_CONTENT, response.status_code)
        self.assertFalse(Plate.objects.filter(barcode="EXP_1").exists())
        # Only the transfer into EXP_2 is left
        self.assertEqual(
            [30], list(WellWithdrawal.objects.values_list("amount", flat=True))
        )
        self.assertEqual(970, self.library_well.amount)

    def test_a_plate_deleted_and_mapped_again_takes_the_volume_once(self):
        plate = self.mapped_plate("EXP_1")

        self.client.delete(f"/api/plates/{plate.id}/")
        self.mapped_plate("EXP_1")

        self.assertEqual(980, self.library_well.amount)

    def test_deleting_a_plate_in_the_admin_deletes_the_withdrawals_into_it(self):
        admin = PlateAdmin(Plate, AdminSite())
        first = self.mapped_plate("EXP_1")
        self.mapped_plate("EXP_2")
        self.mapped_plate("EXP_3")

        admin.delete_model(request=None, obj=first)
        admin.delete_queryset(
            request=None, queryset=Plate.objects.filter(barcode="EXP_2")
        )

        self.assertEqual(
            ["EXP_3"],
            list(
                Plate.objects.filter(experiment=self.experiment).values_list(
                    "barcode", flat=True
                )
            ),
        )
        self.assertEqual(1, WellWithdrawal.objects.count())
        self.assertEqual(980, self.library_well.amount)

    def test_deleting_an_experiment_keeps_the_withdrawals(self):
        self.mapped_plate("EXP_1")

        self.experiment.delete()

        withdrawal = WellWithdrawal.objects.get()
        self.assertIsNone(withdrawal.target_well)
        self.assertEqual(980, self.library_well.amount)
