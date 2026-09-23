"""
Tests for the fill_sdf_amounts command: it fills the amounts of an imported
SDF library and changes nothing else.
"""

import shutil
import tempfile
from os.path import join

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase
from rdkit import Chem

from compoundlib.models import Compound, CompoundLibrary
from core.models import Plate, PlateDimension, WellCompound

MAPPING = (
    "compound:\n  identifier: ID\n  name: NAME\n  structure: Structure\n"
    "plate:\n  barcode: [Barcode_Copy1, Barcode_Copy2]\n"
    "  position: POS_IN_PLATE\n  amount: [Vol_Copy1, Vol_Copy2]\n"
)


class FillSdfAmountsTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.folder)
        self.mapping_file = join(self.folder, "mapping.yml")
        with open(self.mapping_file, "w") as file:
            file.write(MAPPING)
        self.sdf_file = self.write_sdf(
            {"NAME": "Aspirin", "POS_IN_PLATE": "A1", "Vol_Copy1": "24.0"},
            {"NAME": "Caffeine", "POS_IN_PLATE": "B2", "Vol_Copy1": "6"},
        )
        # A library as the SDF import stored it before September 2026: all amounts 0
        call_command(
            "import",
            "sdf",
            input_file=self.sdf_file,
            mapping_file=self.mapping_file,
            library_name="Library",
        )
        WellCompound.objects.update(amount=0)

    def write_sdf(self, *records):
        """One molecule per record, every record on the plates COPY_1 and COPY_2."""
        path = join(self.folder, "library.sdf")
        writer = Chem.SDWriter(path)
        for properties in records:
            molecule = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
            properties = {
                "Barcode_Copy1": "COPY_1",
                "Barcode_Copy2": "COPY_2",
                "Vol_Copy2": "<24",
                **properties,
            }
            for key, value in properties.items():
                molecule.SetProp(key, value)
            writer.write(molecule)
        writer.close()
        return path

    def fill(self, **options):
        with self.assertLogs("API", level="INFO") as logs:
            call_command(
                "fill_sdf_amounts",
                library_name="Library",
                input_file=self.sdf_file,
                mapping_file=self.mapping_file,
                **options,
            )
        return "\n".join(logs.output)

    def amounts(self):
        """e.g. [("COPY_1", "A1", "Aspirin", 24000.0), ...]"""
        return [
            (
                well_compound.well.plate.barcode,
                well_compound.well.hr_position,
                well_compound.compound.name,
                well_compound.amount,
            )
            for well_compound in WellCompound.objects.order_by(
                "well__plate__barcode", "well__position"
            )
        ]

    def test_the_amounts_are_filled_in_nanoliter(self):
        compounds_before = list(Compound.objects.values())

        output = self.fill()

        self.assertEqual(
            [
                ("COPY_1", "A1", "Aspirin", 24000.0),
                ("COPY_1", "B2", "Caffeine", 6000.0),
                ("COPY_2", "A1", "Aspirin", 0.0),
                ("COPY_2", "B2", "Caffeine", 0.0),
            ],
            self.amounts(),
        )
        self.assertEqual(compounds_before, list(Compound.objects.values()))
        self.assertIn("Library Library: 2 amounts changed, 2 already", output)
        self.assertIn("Column Vol_Copy2 has no exact volume ('<24' in 2 wells)", output)

    def test_a_second_run_changes_nothing(self):
        self.fill()

        output = self.fill()

        self.assertIn("0 amounts changed, 4 already had the amount", output)

    def test_a_dry_run_stores_nothing(self):
        output = self.fill(dry_run=True)

        self.assertEqual({0.0}, set(WellCompound.objects.values_list("amount", flat=True)))
        self.assertIn("Library Library: 2 amounts changed", output)
        self.assertIn("Dry run: nothing was stored.", output)

    def test_a_plate_of_another_library_is_not_changed(self):
        other_library = CompoundLibrary.objects.create(name="Other")
        Plate.objects.filter(barcode="COPY_1").update(library=other_library)

        output = self.fill()

        self.assertEqual(
            {0.0}, set(WellCompound.objects.values_list("amount", flat=True))
        )
        self.assertEqual(
            "Other", Plate.objects.get(barcode="COPY_1").library.name
        )
        self.assertIn(
            "1 plates of the file are not in library Library and were skipped: COPY_1.",
            output,
        )

    def test_a_column_that_is_not_a_volume_is_not_used(self):
        self.sdf_file = self.write_sdf(
            {"NAME": "Aspirin", "POS_IN_PLATE": "A1", "Vol_Copy1": "24.0",
             "Weight_Copy2": "247.0"},
            {"NAME": "Caffeine", "POS_IN_PLATE": "B2", "Vol_Copy1": "6",
             "Weight_Copy2": "262.0"},
        )
        with open(self.mapping_file, "w") as file:
            file.write(MAPPING.replace("Vol_Copy2", "Weight_Copy2"))

        output = self.fill()

        self.assertEqual(
            [0.0, 0.0],
            list(
                WellCompound.objects.filter(well__plate__barcode="COPY_2")
                .values_list("amount", flat=True)
            ),
        )
        self.assertIn("2 amounts changed", output)
        self.assertIn(
            "The amounts in column Weight_Copy2 are not stored (set to 0)", output
        )

    def test_a_plate_copy_without_barcode_is_skipped(self):
        self.sdf_file = self.write_sdf(
            {"NAME": "Aspirin", "POS_IN_PLATE": "A1", "Vol_Copy1": "24.0",
             "Barcode_Copy2": ""},
        )
        # Like plate 592 of ActiTarg: a plate with an empty barcode in the library
        Plate.objects.create(
            barcode="",
            dimension=PlateDimension.objects.get(name="dim_384_16x24"),
            library=CompoundLibrary.objects.get(name="Library"),
        )

        output = self.fill()

        self.assertIn(
            "1 plate copies of the file have no barcode and were skipped.", output
        )
        self.assertIn("(no barcode)", output)

    def test_an_unknown_library_is_an_error(self):
        with self.assertRaisesMessage(CommandError, "There is no library Missing."):
            call_command(
                "fill_sdf_amounts",
                library_name="Missing",
                input_file=self.sdf_file,
                mapping_file=self.mapping_file,
            )
