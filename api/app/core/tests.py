from django.test import TestCase

from compoundlib.models import Compound
from .models import Plate, PlateDimension, Well, WellWithdrawal, WellCompound
from .mapping import Mapping, MappingList
from .helper import charToAlphaPos


class MappingTest(TestCase):
    def test_charToAlphaPos(self):
        self.assertEqual(1, charToAlphaPos('A'))
        self.assertEqual(16, charToAlphaPos('P'))
        self.assertEqual(17, charToAlphaPos('q'))
        self.assertEqual(26, charToAlphaPos('z'))
        self.assertRaises(Exception, lambda: charToAlphaPos('bla'))

    def test_positionMapping(self):
        plateDimension = PlateDimension.objects.create(name='bla', rows=2, cols=3)
        pos = plateDimension.position('B2')
        self.assertEqual(4, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual('B2', hr_pos)

        plateDimension = PlateDimension.objects.create(name='bla', rows=4, cols=5)
        pos = plateDimension.position('C3')
        self.assertEqual(12, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual('C3', hr_pos)

        pos = plateDimension.position('D2')
        self.assertEqual(16, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual('D2', hr_pos)

        plateDimension = PlateDimension.objects.create(name='bla', rows=16, cols=22)
        pos = plateDimension.position('A03')
        self.assertEqual(2, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual('A3', hr_pos)


class PlateTest(TestCase):
    def setUp(self):
        self.dimension = PlateDimension.objects.create(
            name='dim_3x2',
            cols=3,
            rows=2
        )
        self.sourcePlate = Plate.objects.create(
            barcode='000001',
            dimension=self.dimension,
        )
        for i in range(self.dimension.num_wells):
            comp = Compound.objects.create(
                identifier=f"comp{i}",
                structure='c=c'
            )
            well = Well.objects.create(
                plate=self.sourcePlate,
                position=i
            )
            WellCompound.objects.create(
                well=well,
                compound=comp,
                amount=1
            )
        self.targetPlate = Plate.objects.create(
            barcode='000002',
            dimension=self.dimension
        )

    def test_plate_mapping(self):
        """Mapping a plate to another"""
        mappingList = MappingList()
        mappingList.add(Mapping(0, 5))
        mappingList.add(Mapping(1, 4))
        mappingList.add(Mapping(2, 3))
        mappingList.add(Mapping(3, 2))
        mappingList.add(Mapping(4, 1))
        mappingList.add(Mapping(5, 0))
        self.sourcePlate.map(mappingList, self.targetPlate)

        targetWells = self.targetPlate.wells.all().order_by('position')
        self.assertEqual('comp5', targetWells[0].compounds.first().identifier)
        self.assertEqual('comp3', targetWells[2].compounds.first().identifier)
        self.assertEqual(1, len(targetWells[3].donors.all()))
        self.assertEqual(2, targetWells[3].donors.first().well.position)
        self.assertEqual('comp2', targetWells[3].donors.first().well.compounds.first().identifier)

    def test_plate_mapping_with_amounts(self):
        """Mapping a plate to another"""
        mappingList = MappingList()
        mappingList.add(Mapping(0, 5, 0.1))
        mappingList.add(Mapping(1, 4, 0.1))
        mappingList.add(Mapping(2, 3, 0.1))
        mappingList.add(Mapping(3, 2, 0.1))
        mappingList.add(Mapping(4, 1, 0.1))
        mappingList.add(Mapping(5, 0, 0.1))
        self.sourcePlate.map(mappingList, self.targetPlate)

        sourceWells = self.sourcePlate.wells.all().order_by('position')
        for sourceWell in sourceWells:
            self.assertEqual(0.9, sourceWell.amount)

        targetWells = self.targetPlate.wells.all().order_by('position')
        for targetWell in targetWells:
            self.assertEqual(0.1, targetWell.amount)

    def test_plate_copy(self):
        """Copy a plate to another"""

        mappingList = MappingList.one_to_one(self.dimension.num_wells)
        self.sourcePlate.map(mappingList, self.targetPlate)

        targetWells = self.targetPlate.wells.all().order_by('position')
        for p in range(len(targetWells)):
            self.assertEqual(f"comp{p}", targetWells[p].compounds.first().identifier)


class WellTest(TestCase):
    def setUp(self):
        self.dimension = PlateDimension.objects.create(
            name='dim_3x2',
            cols=3,
            rows=2
        )
        self.plate = Plate.objects.create(
            barcode='123456789',
            dimension=self.dimension,
        )

    def test_withdrawal(self):
        compound = Compound.objects.create(identifier='ABC', structure='A-B-C')
        well = Well.objects.create(position=0, plate=self.plate)
        WellCompound.objects.create(well=well, compound=compound, amount=100)
        WellWithdrawal.objects.create(well=well, amount=10)
        WellWithdrawal.objects.create(well=well, amount=10)
        self.assertEqual(80, well.amount)

    def test_source_plate_discovery(self):
        """ Test if we can find the correct source plate if we map twice """

        def __fillPlateWithCompounds(plate):
            for i in range(plate.dimension.num_wells):
                comp = Compound.objects.create(
                    identifier=f"{plate.barcode}_comp{i}",
                    structure=f"comp{i}"
                )
                well = Well.objects.create(
                    plate=plate,
                    position=i
                )
                WellCompound.objects.create(
                    well=well,
                    compound=comp,
                    amount=1
                )

        plate1 = Plate.objects.create(
            barcode='0001',
            dimension=self.dimension,
        )

        plate2 = Plate.objects.create(
            barcode='0002',
            dimension=self.dimension,
        )

        plate3 = Plate.objects.create(
            barcode='0003',
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
