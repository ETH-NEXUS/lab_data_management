import csv
from django.db import IntegrityError, transaction
from django.test import TestCase
from compoundlib.models import Compound, CompoundLibrary
from core.models import (
    Plate,
    PlateDimension,
    Well,
    WellCompound,
    Experiment,
    Project,
)
from platetemplate.models import PlateTemplate, PlateTemplateCategory
from core.mapping import Mapping, MappingList


class PlateTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        self.sourcePlate = Plate.objects.create(
            barcode="000001",
            dimension=self.dimension,
        )
        for i in range(self.dimension.num_wells):
            comp = Compound.objects.create(name=f"comp{i}", structure="c=c")
            well = Well.objects.create(plate=self.sourcePlate, position=i)
            WellCompound.objects.create(well=well, compound=comp, amount=1)
        self.targetPlate = Plate.objects.create(
            barcode="000002", dimension=self.dimension
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

        targetWells = self.targetPlate.wells.all().order_by("position")
        self.assertEqual("comp5", targetWells[0].compounds.first().name)
        self.assertEqual("comp3", targetWells[2].compounds.first().name)
        self.assertEqual(1, len(targetWells[3].donors.all()))
        self.assertEqual(2, targetWells[3].donors.first().well.position)
        self.assertEqual(
            "comp2",
            targetWells[3].donors.first().well.compounds.first().name,
        )

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

        sourceWells = self.sourcePlate.wells.all().order_by("position")
        for sourceWell in sourceWells:
            self.assertEqual(0.9, sourceWell.amount)

        targetWells = self.targetPlate.wells.all().order_by("position")
        for targetWell in targetWells:
            self.assertEqual(0.1, targetWell.amount)

    def test_plate_copy(self):
        """Copy a plate to another"""

        mappingList = MappingList.one_to_one(self.dimension.num_wells)
        self.sourcePlate.map(mappingList, self.targetPlate)

        targetWells = self.targetPlate.wells.all().order_by("position")
        for p in range(len(targetWells)):
            self.assertEqual(f"comp{p}", targetWells[p].compounds.first().name)

    def test_plate_mapping_from_csv(self):
        """CSV Mapping test"""
        delimiter = ";"
        csv_file = "test.csv"
        with open(csv_file, "w", newline="") as cf:
            writer = csv.DictWriter(cf, ("from", "to", "amount"), delimiter=delimiter)
            writer.writeheader()
            writer.writerow({"from": 0, "to": 5, "amount": 10})
            writer.writerow({"from": 1, "to": 4, "amount": 20})

        mappingList = MappingList.from_csv(
            csv_file, "from", "to", "amount", delimiter=";"
        )
        self.assertEqual(0, mappingList[0].from_pos)
        self.assertEqual(5, mappingList[0].to_pos)
        self.assertEqual(1, mappingList[1].from_pos)
        self.assertEqual(4, mappingList[1].to_pos)

        self.sourcePlate.map(mappingList, self.targetPlate)
        self.assertEqual(
            self.sourcePlate.wells.get(position=0).compounds.first(),
            self.targetPlate.wells.get(position=5).compounds.first(),
        )

    def test_plate_mapping_from_csv_with_illegal_values(self):
        """CSV Mapping test with illegal values"""
        delimiter = ";"
        csv_file = "test.csv"
        with open(csv_file, "w", newline="") as cf:
            writer = csv.DictWriter(cf, ("from", "to", "amount"), delimiter=delimiter)
            writer.writeheader()
            writer.writerow({"from": "a", "to": "c", "amount": 10})
            writer.writerow({"from": "b", "to": "d", "amount": 20})

        self.assertRaises(
            ValueError,
            MappingList.from_csv,
            csv_file,
            "from",
            "to",
            "amount",
            delimiter=";",
        )

    def test_only_library_or_experiment_or_template_constraint(self):
        """Check the constraint that only library or experiment can have a value on a plate"""
        library = CompoundLibrary.objects.create(name="CL")
        project = Project.objects.create(name="Proj")
        experiment = Experiment.objects.create(name="Exp", project=project)
        category = PlateTemplateCategory.objects.create(name="PTC")
        template = PlateTemplate.objects.create(name="PT", category=category)
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode="bar1",
                    dimension=self.dimension,
                    library=library,
                    experiment=experiment,
                    template=template,
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode="bar2",
                    dimension=self.dimension,
                    library=library,
                    experiment=experiment,
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode="bar3",
                    dimension=self.dimension,
                    experiment=experiment,
                    template=template,
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass
        try:
            with transaction.atomic():
                Plate.objects.create(
                    barcode="bar4",
                    dimension=self.dimension,
                    library=library,
                    template=template,
                )
            self.fail("Should raise an IntegrityError")
        except IntegrityError:
            pass

    def test_create_library_plate(self):
        """Check the creation of a library plate"""
        library = CompoundLibrary.objects.create(name="CL")
        plate = Plate.objects.create(
            barcode="123456", dimension=self.dimension, library=library
        )
        self.assertIsNotNone(plate)

    def test_create_experiment_plate(self):
        """Check the creation of an experiment plate"""
        project = Project.objects.create(name="Proj")
        experiment = Experiment.objects.create(name="Exp", project=project)
        plate = Plate.objects.create(
            barcode="123456", dimension=self.dimension, experiment=experiment
        )
        self.assertIsNotNone(plate)
