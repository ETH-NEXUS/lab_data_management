import csv

from django.db import IntegrityError, transaction
from django.test import TestCase

from compoundlib.models import Compound, CompoundLibrary
from .models import Plate, PlateDimension, Well, WellWithdrawal, WellCompound, Experiment, Project, MeasurementFeature, Measurement, WellType
from platetemplate.models import PlateTemplate, PlateTemplateCategory
from .helper import charToAlphaPos
from .mapping import Mapping, MappingList
from .helper import charToAlphaPos

import csv
import numpy as np


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

    def test_plate_mapping_from_csv(self):
        """CSV Mapping test"""
        delimiter = ';'
        csv_file = 'test.csv'
        with open(csv_file, 'w', newline='') as cf:
            writer = csv.DictWriter(cf, ('from', 'to', 'amount'), delimiter=delimiter)
            writer.writeheader()
            writer.writerow({'from': 0, 'to': 5, 'amount': 10})
            writer.writerow({'from': 1, 'to': 4, 'amount': 20})

        mappingList = MappingList.from_csv(csv_file, 'from', 'to', 'amount', delimiter=';')
        self.assertEqual(0, mappingList[0].from_pos)
        self.assertEqual(5, mappingList[0].to_pos)
        self.assertEqual(1, mappingList[1].from_pos)
        self.assertEqual(4, mappingList[1].to_pos)

        self.sourcePlate.map(mappingList, self.targetPlate)
        self.assertEqual(
            self.sourcePlate.wells.get(
                position=0).compounds.first(), self.targetPlate.wells.get(
                position=5).compounds.first())

    def test_plate_mapping_from_csv_with_illegal_values(self):
        """CSV Mapping test with illegal values"""
        delimiter = ';'
        csv_file = 'test.csv'
        with open(csv_file, 'w', newline='') as cf:
            writer = csv.DictWriter(cf, ('from', 'to', 'amount'), delimiter=delimiter)
            writer.writeheader()
            writer.writerow({'from': 'a', 'to': 'c', 'amount': 10})
            writer.writerow({'from': 'b', 'to': 'd', 'amount': 20})

        self.assertRaises(ValueError, MappingList.from_csv, csv_file, 'from', 'to', 'amount', delimiter=';')

    def test_only_library_or_experiment_or_template_constraint(self):
        """ Check the constraint that only library or experiment can have a value on a plate"""
        library = CompoundLibrary.objects.create(name='CL')
        project = Project.objects.create(name='Proj')
        experiment = Experiment.objects.create(name='Exp', project=project)
        category = PlateTemplateCategory.objects.create(name='PTC')
        template = PlateTemplate.objects.create(name='PT', category=category)
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode='bar1',
                    dimension=self.dimension,
                    library=library,
                    experiment=experiment,
                    template=template
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode='bar2',
                    dimension=self.dimension,
                    library=library,
                    experiment=experiment
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode='bar3',
                    dimension=self.dimension,
                    experiment=experiment,
                    template=template
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode='bar4',
                    dimension=self.dimension,
                    library=library,
                    template=template
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass

    def test_create_library_plate(self):
        """ Check the creation of a library plate"""
        library = CompoundLibrary.objects.create(name='CL')
        plate = Plate.objects.create(
            barcode='123456',
            dimension=self.dimension,
            library=library
        )
        self.assertIsNotNone(plate)

    def test_create_experiment_plate(self):
        """ Check the creation of an experiment plate"""
        project = Project.objects.create(name='Proj')
        experiment = Experiment.objects.create(name='Exp', project=project)
        plate = Plate.objects.create(
            barcode='123456',
            dimension=self.dimension,
            experiment=experiment
        )
        self.assertIsNotNone(plate)


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


class StatisticsTest(TestCase):
    fixtures = ('well_types',)

    def setUp(self):
        """
        Creates a plate 4x3 with the following types:values.

        P:0 C:1 C:2 N:3
        P:4 C:5 C:6 N:7
        P:8 C:9 C:10 N:11

        """
        self.dimension = PlateDimension.objects.create(
            name='dim_4x3',
            cols=4,
            rows=3
        )
        self.plate = Plate.objects.create(
            barcode='000001',
            dimension=self.dimension,
        )
        self.measurement_feature = MeasurementFeature.objects.create(
            name='TEST',
            abbrev='TST',
            unit='t'
        )
        for i in range(self.dimension.num_wells):
            comp = Compound.objects.create(
                identifier=f"comp{i}",
                structure='c=c'
            )
            well = Well.objects.create(
                plate=self.plate,
                position=i,
                type=WellType.by_name('P') if i in [0, 4, 8] else WellType.by_name('N') if i in [3, 7, 11] else WellType.by_name('C')
            )
            WellCompound.objects.create(
                well=well,
                compound=comp,
                amount=1
            )
            Measurement.objects.create(
                well=well,
                feature=self.measurement_feature,
                value=i
            )

    def test_z_prime_calculation(self):
        """ Test calculation of z' (prime) factor"""
        expected_z_prime = 1 - (3 * (np.std([0, 4, 8]) + np.std([3, 7, 11])) / abs(np.mean([0, 4, 8]) - np.mean([3, 7, 11])))
        self.assertEqual(expected_z_prime, self.plate.z_prime('TST'))

    def test_z_factor_calculation(self):
        """ Test calculation of z factor """
        expected_z_prime = 1 - (3 * (np.std([0, 4, 8]) + np.std([1, 2, 5, 6, 9, 10])) / abs(np.mean([0, 4, 8]) - np.mean([1, 2, 5, 6, 9, 10])))
        self.assertEqual(expected_z_prime, self.plate.z_factor('TST'))

    def test_z_score_calculation(self):
        """ Test the calculation of the z score per well """
        c_pos = [1, 2, 5, 6, 9, 10]
        mean = np.mean(c_pos)
        std = np.std(c_pos)
        expected_scores = [(v - mean) / std for v in range(self.plate.num_wells)]
        calculated_scores = list(filter(lambda v: v is not None, self.plate.z_scores('TST')))
        self.assertEqual(len(expected_scores), len(calculated_scores))
        self.assertListEqual(expected_scores, calculated_scores)
