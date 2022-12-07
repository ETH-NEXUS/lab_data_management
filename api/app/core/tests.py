from django.test import TestCase

from compoundlib.models import Compound
from .models import Plate, PlateDimension, Well
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
                structure='c=c',
                smile='c=c'
            )
            Well.objects.create(
                plate=self.sourcePlate,
                position=i,
                compound=comp
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
        self.sourcePlate.map(mappingList, self.targetPlate, 1.0)

        targetWells = self.targetPlate.wells.all().order_by('position')
        self.assertEqual('comp5', targetWells[0].compound.identifier)
        self.assertEqual('comp3', targetWells[2].compound.identifier)
        self.assertEqual(1, len(targetWells[3].from_well.all()))
        self.assertEqual(2, targetWells[3].from_well.first().position)
        self.assertEqual('comp2', targetWells[3].from_well.first().compound.identifier)

    def test_plate_copy(self):
        """Copy a plate to another"""

        mappingList = MappingList.one_to_one(self.dimension.num_wells)
        self.sourcePlate.map(mappingList, self.targetPlate, 1.0)

        targetWells = self.targetPlate.wells.all().order_by('position')
        for p in range(len(targetWells)):
            self.assertEqual(f"comp{p}", targetWells[p].compound.identifier)
