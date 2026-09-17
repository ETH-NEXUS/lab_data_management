"""
Tests for the errors of the import command, as the management page shows them.
"""

import shutil
import tempfile
from os.path import join
from unittest import mock

from django.core.cache import cache
from django.test import TestCase
from django.urls import reverse
from rdkit import Chem

from compoundlib.models import Compound
from importer.mapping import SdfMapping
from core.models import Plate, Project, WellCompound

# A library plate file: the compound names, one empty line, the well types
LIBRARY_PLATE_CSV = (
    "Aspirin,null,Caffeine\nnull,null,null\n\nC,null,P\nnull,null,null\n"
)


class ImportCommandTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        cache.clear()
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def write(self, name, text):
        path = join(self.folder, name)
        with open(path, "w") as file:
            file.write(text)
        return path

    def write_sdf(self, properties):
        """One molecule with the given SDF properties, e.g. {"NAME": "Aspirin"}."""
        path = join(self.folder, "library.sdf")
        molecule = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        for key, value in properties.items():
            molecule.SetProp(key, value)
        writer = Chem.SDWriter(path)
        writer.write(molecule)
        writer.close()
        return path

    def run_import(self, what, **form_data):
        data = {"command": "import", "what": what, "room_name": "room_1"}
        data.update(form_data)
        self.client.post(
            reverse("run_command"), {"form_data": data}, content_type="application/json"
        )
        url = reverse("long_polling", args=["room_1"])
        return self.client.get(f"{url}?since=0").json()

    def assertFailedWith(self, output, text):
        """The command failed, and `text` is its only error besides "Command failed."."""
        self.assertEqual("failed", output["status"])
        errors = [
            message["text"]
            for message in output["messages"]
            if message["level"] == "error"
        ]
        self.assertEqual([text, "Command failed."], errors)

    def test_a_library_plate_is_imported(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertEqual("completed", output["status"])
        plate = Plate.objects.get(barcode="LIB_1")
        self.assertEqual(
            [(0, "Aspirin", "C"), (2, "Caffeine", "P")],
            list(
                WellCompound.objects.filter(well__plate=plate)
                .order_by("well__position")
                .values_list("well__position", "compound__name", "well__type__name")
            ),
        )

    def test_a_library_plate_file_saved_by_excel_uses_the_existing_compound(self):
        # Excel writes an invisible BOM character before the first compound name
        aspirin = Compound.objects.create(name="Aspirin")
        path = self.write("plate.csv", "\ufeff" + LIBRARY_PLATE_CSV)

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertEqual("completed", output["status"])
        self.assertEqual(
            ["Aspirin", "Caffeine"],
            sorted(Compound.objects.values_list("name", flat=True)),
        )
        self.assertEqual(
            aspirin,
            WellCompound.objects.get(
                well__plate__barcode="LIB_1", well__position=0
            ).compound,
        )

    def test_a_library_plate_file_that_does_not_exist(self):
        output = self.run_import(
            "library_plate",
            input_file="/no/such/plate.csv",
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertFailedWith(output, "File does not exist: /no/such/plate.csv")
        self.assertFalse(Plate.objects.exists())

    def test_a_library_plate_file_with_a_wrong_format_says_why(self):
        path = self.write("plate.csv", "Aspirin,null\nnull,null\n")

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertFailedWith(
            output,
            "File format is incorrect: The file should contain exactly one empty line.",
        )

    def test_a_control_plate_for_a_project_that_does_not_exist(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)

        output = self.run_import(
            "library_plate",
            input_file=path,
            project_name="No such project",
            plate_barcode="CTRL_1",
            is_control_plate=True,
        )

        self.assertFailedWith(output, "Project No such project does not exist.")
        self.assertFalse(Plate.objects.exists())

    def test_a_library_plate_without_library_and_project(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)

        output = self.run_import(
            "library_plate", input_file=path, plate_barcode="LIB_1"
        )

        self.assertFailedWith(
            output, "Please specify a library name or a project name."
        )

    def test_a_plate_barcode_that_already_exists(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)
        Project.objects.create(name="Project")
        self.run_import(
            "library_plate", input_file=path, project_name="Project", plate_barcode="P1"
        )

        output = self.run_import(
            "library_plate", input_file=path, project_name="Project", plate_barcode="P1"
        )

        self.assertFailedWith(output, "Plate with barcode P1 already exists.")

    def test_an_sdf_library_is_imported(self):
        path = self.write_sdf(
            {
                "NAME": "Aspirin",
                "PLATE_NUMBER1": "SDF_1",
                "POS_IN_PLATE": "A1",
                "PLATE_AMOUNT1": "10",
            }
        )

        output = self.run_import("sdf", input_file=path, library_name="Library")

        self.assertEqual("completed", output["status"])
        well_compound = WellCompound.objects.get()
        self.assertEqual(
            ("SDF_1", 0, "Aspirin"),
            (
                well_compound.well.plate.barcode,
                well_compound.well.position,
                well_compound.compound.name,
            ),
        )

    def test_an_sdf_file_that_does_not_exist(self):
        output = self.run_import("sdf", input_file="/no/such/library.sdf")

        self.assertFailedWith(output, "File does not exist: /no/such/library.sdf")

    def test_an_sdf_file_without_the_mapped_columns(self):
        path = self.write_sdf({"NAME": "Aspirin", "POS_IN_PLATE": "A1"})

        output = self.run_import("sdf", input_file=path, library_name="Library")

        self.assertFailedWith(
            output,
            f"These columns are not in the SDF file {path}: PLATE_NUMBER1, "
            "PLATE_AMOUNT1. Check the mapping file.",
        )
        self.assertFalse(Compound.objects.exists())

    def test_a_mapping_file_with_a_wrong_structure(self):
        path = self.write_sdf({"NAME": "Aspirin"})
        mapping_file = self.write("mapping.yml", "compound:\n  name: NAME\n")

        output = self.run_import("sdf", input_file=path, mapping_file=mapping_file)

        self.assertEqual("failed", output["status"])
        self.assertTrue(
            output["messages"][0]["text"].startswith(
                "MappingFileSchemaError: Error in mapping file schema."
            )
        )

    def test_an_unexpected_error_shows_its_type_and_is_logged_with_traceback(self):
        with mock.patch(
            "importer.management.commands.import.Command.sdf",
            side_effect=KeyError("NAME"),
        ):
            with self.assertLogs("API", level="ERROR") as logs:
                output = self.run_import("sdf", input_file="/data/library.sdf")

        self.assertFailedWith(output, "KeyError: 'NAME'")
        failure = logs.records[-1]
        self.assertEqual("Command import sdf failed", failure.getMessage())
        self.assertIsNotNone(failure.exc_info)

    def test_a_mapping_file_does_not_change_the_default_mapping(self):
        mapping_file = self.write(
            "mapping.yml",
            "compound:\n  identifier: ID\n  name: CompoundName\n  structure: MOL\n"
            "plate:\n  barcode: Barcode\n  position: Well\n  amount: Volume\n",
        )

        self.assertEqual("CompoundName", SdfMapping(mapping_file).name)
        self.assertEqual("NAME", SdfMapping().name)
