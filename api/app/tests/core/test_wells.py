from django.test import TestCase
from compoundlib.models import Compound
from core.models import (
    Plate,
    PlateDimension,
    Well,
    WellWithdrawal,
    WellCompound,
)


class WellTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        self.plate = Plate.objects.create(
            barcode="123456789",
            dimension=self.dimension,
        )

    def test_withdrawal(self):
        compound = Compound.objects.create(name="ABC", structure="A-B-C")
        well = Well.objects.create(position=0, plate=self.plate)
        WellCompound.objects.create(well=well, compound=compound, amount=100)
        WellWithdrawal.objects.create(well=well, amount=10)
        WellWithdrawal.objects.create(well=well, amount=10)
        self.assertEqual(80, well.amount)

    def test_source_plate_discovery(self):
        """Test if we can find the correct source plate if we map twice"""

        def __fillPlateWithCompounds(plate):
            for i in range(plate.dimension.num_wells):
                comp = Compound.objects.create(
                    name=f"{plate.barcode}_comp{i}", structure=f"comp{i}"
                )
                well = Well.objects.create(plate=plate, position=i)
                WellCompound.objects.create(well=well, compound=comp, amount=1)

        plate1 = Plate.objects.create(
            barcode="0001",
            dimension=self.dimension,
        )

        plate2 = Plate.objects.create(
            barcode="0002",
            dimension=self.dimension,
        )

        plate3 = Plate.objects.create(
            barcode="0003",
            dimension=self.dimension,
        )

        __fillPlateWithCompounds(plate1)

        plate1.copy(plate2, 0.6)

        for well in plate1.wells.all():
            self.assertEqual(0.4, well.amount)

        for well in plate2.wells.all():
            self.assertEqual(0.6, well.amount)
            self.assertEqual(plate1, well.donors.first().well.plate)

        plate2.copy(plate3, 0.2)

        for well in plate2.wells.all():
            self.assertEqual(0.4, well.amount)
            self.assertEqual(plate1, well.donors.first().well.plate)

        for well in plate3.wells.all():
            self.assertEqual(0.2, well.amount)
            self.assertEqual(plate2, well.donors.first().well.plate)
            self.assertEqual(plate1, well.donors.first().well.donors.first().well.plate)
