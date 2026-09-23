"""
Tests for the errors of the import command, as the management page shows them.
"""

from importlib import import_module
from os.path import join
from unittest import mock

from django.test import TestCase
from rdkit import Chem

from compoundlib.models import Compound, CompoundLibrary
from tests.management.base import ManagementPageTestCase
from importer.mapping import SdfMapping
from core.models import Plate, Project, Well, WellCompound

# The command module is called "import", a Python keyword, so it cannot be imported directly
yes_or_no = import_module("importer.management.commands.import").yes_or_no

# A library plate file: the compound names, one empty line, the well types
LIBRARY_PLATE_CSV = (
    "Aspirin,null,Caffeine\nnull,null,null\n\nC,null,P\nnull,null,null\n"
)


class ImportCommandTest(ManagementPageTestCase):
    def write_sdf(self, *records):
        """One molecule per record of SDF properties, e.g. {"NAME": "Aspirin"}."""
        path = join(self.folder, "library.sdf")
        writer = Chem.SDWriter(path)
        for properties in records:
            molecule = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
            for key, value in properties.items():
                molecule.SetProp(key, value)
            writer.write(molecule)
        writer.close()
        return path

    def run_import(self, what, **form_data):
        """Starts an import command and returns its output."""
        self.start_command(command="import", what=what, **form_data)
        return self.read_output()

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

    def test_an_unknown_kind_of_import_is_an_error(self):
        output = self.run_import("plates", input_file=join(self.folder, "plate.csv"))

        self.assertEqual("failed", output["status"])
        self.assertIn("invalid choice", self.errors(output)[0])

    def test_a_library_plate_is_not_a_control_plate(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)

        self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
            is_control_plate=False,
        )

        self.assertIs(False, Plate.objects.get(barcode="LIB_1").is_control_plate)

    def test_a_control_plate_of_a_project(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)
        project = Project.objects.create(name="Project")

        self.run_import(
            "library_plate",
            input_file=path,
            project_name="Project",
            plate_barcode="CTRL_1",
            is_control_plate=True,
        )

        plate = Plate.objects.get(barcode="CTRL_1")
        self.assertIs(True, plate.is_control_plate)
        self.assertEqual(project, plate.project)

    def test_an_error_in_the_middle_of_a_plate_file_stores_nothing(self):
        # Aspirin is imported first, then the well type "XX" does not exist
        path = self.write(
            "plate.csv", LIBRARY_PLATE_CSV.replace("C,null,P", "C,null,XX")
        )

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertEqual("failed", output["status"])
        self.assertTrue(
            self.errors(output)[0].startswith(
                "Unknown well type 'XX' in well A3. Known well types: C, P, N, R1,"
            )
        )
        self.assertIn(
            {
                "level": "warning",
                "text": "This import failed, nothing of it was stored. The reason:",
            },
            output["messages"],
        )
        self.assertFalse(CompoundLibrary.objects.exists())
        self.assertFalse(Plate.objects.exists())
        self.assertFalse(Well.objects.exists())
        self.assertFalse(Compound.objects.exists())

    def test_an_error_in_the_middle_of_an_sdf_file_stores_nothing(self):
        record = {"NAME": "Aspirin", "PLATE_NUMBER1": "SDF_1", "PLATE_AMOUNT1": "10"}
        path = self.write_sdf(
            {**record, "POS_IN_PLATE": "A1"},
            {**record, "NAME": "Caffeine", "POS_IN_PLATE": "not a well"},
        )

        output = self.run_import("sdf", input_file=path, library_name="Library")

        self.assertEqual("failed", output["status"])
        self.assertFalse(CompoundLibrary.objects.exists())
        self.assertFalse(Plate.objects.exists())
        self.assertFalse(Compound.objects.exists())

    def test_a_library_plate_file_that_does_not_exist(self):
        output = self.run_import(
            "library_plate",
            input_file=join(self.folder, "missing.csv"),
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertFailedWith(
            output, f"File does not exist: {join(self.folder, 'missing.csv')}"
        )
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

    def test_a_plate_file_with_a_short_line_is_refused(self):
        path = self.write(
            "plate.csv", LIBRARY_PLATE_CSV.replace("null,null,null", "null,null")
        )

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        self.assertFailedWith(
            output,
            "File format is incorrect: Every line must have 3 cells, like the "
            "first line, but line 2 has 2.",
        )

    def test_the_new_compounds_of_a_plate_are_reported_in_one_message(self):
        path = self.write("plate.csv", LIBRARY_PLATE_CSV)

        output = self.run_import(
            "library_plate",
            input_file=path,
            library_name="Library",
            plate_barcode="LIB_1",
        )

        texts = [message["text"] for message in output["messages"]]
        self.assertIn("Created 2 new compounds: Aspirin, Caffeine", texts)

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

    def test_the_sdf_volumes_are_stored_in_nanoliter(self):
        mapping_file = self.write(
            "mapping.yml",
            "compound:\n  identifier: ID\n  name: NAME\n  structure: Structure\n"
            "plate:\n  barcode: [Barcode_Copy1, Barcode_Copy2]\n"
            "  position: POS_IN_PLATE\n  amount: [Vol_Copy1, Vol_Copy2]\n",
        )
        path = self.write_sdf(
            {
                "NAME": "Aspirin",
                "POS_IN_PLATE": "A1",
                "Barcode_Copy1": "COPY_1",
                "Vol_Copy1": "24.0",
                "Barcode_Copy2": "COPY_2",
                "Vol_Copy2": "<24",
            }
        )

        output = self.run_import(
            "sdf", input_file=path, library_name="Library", mapping_file=mapping_file
        )

        self.assertEqual("completed", output["status"])
        self.assertEqual(
            [("COPY_1", 24000.0), ("COPY_2", 0.0)],
            list(
                WellCompound.objects.order_by("well__plate__barcode").values_list(
                    "well__plate__barcode", "amount"
                )
            ),
        )
        self.assertIn(
            {
                "level": "warning",
                "text": "Column Vol_Copy2 has no exact volume ('<24' in 1 wells), "
                "so these amounts are set to 0.",
            },
            output["messages"],
        )

    def test_an_sdf_amount_of_an_unknown_unit_is_not_stored(self):
        path = self.write_sdf(
            {
                "NAME": "Aspirin",
                "PLATE_NUMBER1": "SDF_1",
                "POS_IN_PLATE": "A1",
                "PLATE_AMOUNT1": "247.0",
            }
        )

        output = self.run_import("sdf", input_file=path, library_name="Library")

        self.assertEqual("completed", output["status"])
        self.assertEqual(0, WellCompound.objects.get().amount)
        self.assertIn(
            {
                "level": "warning",
                "text": "The amounts in column PLATE_AMOUNT1 are not stored (set to 0), "
                "because only the Vol_Copy… columns are known to be volumes in µL.",
            },
            output["messages"],
        )

    def test_an_sdf_file_that_does_not_exist(self):
        path = join(self.folder, "missing.sdf")

        output = self.run_import("sdf", input_file=path)

        self.assertFailedWith(output, f"File does not exist: {path}")

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
            self.errors(output)[0].startswith(
                "MappingFileSchemaError: Error in mapping file schema."
            )
        )

    def test_an_unexpected_error_shows_its_type_and_is_logged_with_traceback(self):
        with mock.patch(
            "importer.management.commands.import.Command.sdf",
            side_effect=KeyError("NAME"),
        ):
            with self.assertLogs("API", level="ERROR") as logs:
                output = self.run_import(
                    "sdf", input_file=join(self.folder, "library.sdf")
                )

        self.assertFailedWith(output, "KeyError: 'NAME'")
        failure = logs.records[-1]
        self.assertIn("Command failed:", failure.getMessage())
        self.assertIn("'what': 'sdf'", failure.getMessage())
        self.assertIsNotNone(failure.exc_info)

    def test_a_mapping_file_does_not_change_the_default_mapping(self):
        mapping_file = self.write(
            "mapping.yml",
            "compound:\n  identifier: ID\n  name: CompoundName\n  structure: MOL\n"
            "plate:\n  barcode: Barcode\n  position: Well\n  amount: Volume\n",
        )

        self.assertEqual("CompoundName", SdfMapping(mapping_file).name)
        self.assertEqual("NAME", SdfMapping().name)

    def import_template(self, text):
        path = self.write("template.tsv", text)
        return self.run_import("template", input_file=path, template_name="Screen")

    def test_a_template_keeps_the_whole_well_type_names(self):
        output = self.import_template("C\tR10\tP1\nN\tP\tn4\n")

        self.assertEqual("completed", output["status"])
        plate = Plate.objects.get(barcode="__TEMPL__Default_Screen")
        self.assertEqual((2, 3), (plate.dimension.rows, plate.dimension.cols))
        self.assertEqual(
            ["C", "R10", "P1", "N", "P", "N4"],
            list(plate.wells.order_by("position").values_list("type__name", flat=True)),
        )

    def test_a_template_with_an_unknown_well_type_stores_nothing(self):
        output = self.import_template("C\tC\tC\nC\tR99\tC\n")

        self.assertEqual("failed", output["status"])
        self.assertTrue(
            self.errors(output)[0].startswith(
                "Unknown well type 'R99' in well B2. Known well types:"
            )
        )
        self.assertFalse(Plate.objects.filter(template__isnull=False).exists())

    def test_a_template_with_an_empty_cell_is_refused(self):
        output = self.import_template("C\tC\tC\nC\t\tC\n")

        self.assertFailedWith(output, "Well B2 has no well type in the template file.")

    def test_a_template_with_rows_of_different_length_is_refused(self):
        output = self.import_template("C\tC\tC\nC\tC\n")

        self.assertFailedWith(
            output, "Row 2 of the template file has 2 cells, but row 1 has 3."
        )

    def test_an_empty_sdf_mapping_file_is_a_schema_error(self):
        path = self.write_sdf({"NAME": "Aspirin"})
        mapping_file = self.write("mapping.yml", "")

        output = self.run_import("sdf", input_file=path, mapping_file=mapping_file)

        self.assertTrue(
            self.errors(output)[0].startswith(
                "MappingFileSchemaError: Error in mapping file schema."
            )
        )

    def test_a_missing_template_file_is_an_error_without_refreshing(self):
        path = join(self.folder, "missing.tsv")

        output = self.run_import("template", input_file=path)

        self.assertFailedWith(output, f"File does not exist: {path}")
        texts = [message["text"] for message in output["messages"]]
        self.assertNotIn("Refreshing materialized views...", texts)

    def test_an_sdf_import_into_an_existing_library_says_using(self):
        CompoundLibrary.objects.create(name="Library")
        path = self.write_sdf(
            {
                "NAME": "Aspirin",
                "PLATE_NUMBER1": "SDF_1",
                "POS_IN_PLATE": "A1",
                "PLATE_AMOUNT1": "10",
            }
        )

        output = self.run_import("sdf", input_file=path, library_name="Library")

        texts = [message["text"] for message in output["messages"]]
        self.assertIn("Using library Library.", texts)


class YesOrNoTest(TestCase):
    def test_command_line_values(self):
        self.assertTrue(yes_or_no("yes"))
        self.assertTrue(yes_or_no("True"))
        self.assertFalse(yes_or_no("False"))
        self.assertFalse(yes_or_no("no"))
