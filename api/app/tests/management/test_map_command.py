"""
Tests for the errors of the map command, as the management page shows them.
"""

from os.path import join
from unittest import mock

from core.models import Experiment, Measurement, Project
from importer.mappers import EchoMapper
from tests.importer.test_m1000_parse import ASC_FILE_CONTENT
from tests.management.base import ManagementPageTestCase


class MapCommandTest(ManagementPageTestCase):
    def setUp(self):
        super().setUp()
        project = Project.objects.create(name="Project")
        Experiment.objects.create(name="Experiment", project=project)

    def run_map(self, machine, **form_data):
        """Starts a map command and returns its output."""
        data = {
            "command": "map",
            "machine": machine,
            "path": self.folder,
            "experiment_name": "Experiment",
        }
        data.update(form_data)
        self.start_command(**data)
        return self.read_output()

    def echo_column_file(self, columns):
        lines = [f'{key}: "{name}"' for key, name in columns.items()]
        return self.write("columns.yml", "\n".join(lines) + "\n")

    def test_the_echo_column_file_is_used(self):
        columns = {
            key: f"My {name}" for key, name in EchoMapper.DEFAULT_COLUMNS.items()
        }
        column_file = self.echo_column_file(columns)

        with mock.patch("importer.management.commands.map.EchoMapper.run") as run:
            output = self.run_map("echo", mapping_file=column_file)

        self.assertEqual("completed", output["status"])
        self.assertEqual(columns, run.call_args.kwargs["headers"])

    def test_without_column_file_the_default_echo_columns_are_used(self):
        with mock.patch("importer.management.commands.map.EchoMapper.run") as run:
            self.run_map("echo", mapping_file="")

        self.assertEqual(EchoMapper.DEFAULT_COLUMNS, run.call_args.kwargs["headers"])

    def test_an_echo_column_file_that_does_not_exist(self):
        path = join(self.folder, "missing.yml")

        output = self.run_map("echo", mapping_file=path)

        self.assertFailedWith(output, f"The column file '{path}' could not be found.")

    def test_an_echo_column_file_with_missing_keys(self):
        column_file = self.echo_column_file({"source_well": "Source Well"})

        output = self.run_map("echo", mapping_file=column_file)

        self.assertEqual("failed", output["status"])
        self.assertTrue(
            self.errors(output)[0].startswith(
                f"The column file '{column_file}' is missing these keys: "
                "source_plate_barcode, "
            )
        )

    def test_m1000_without_experiment_name(self):
        output = self.run_map("m1000", experiment_name="")

        self.assertFailedWith(
            output,
            "No experiment name provided. If you would like to add missing plates, "
            "you need to provide the experiment name.",
        )

    def test_an_m1000_value_that_is_not_a_number_stops_before_anything_is_stored(
        self,
    ):
        path = join(self.folder, "20240610-121212_demo_1.asc")
        with open(path, "wb") as file:
            file.write(ASC_FILE_CONTENT.replace(b"A2\tSM1_2\t999", b"A2\tSM1_2\t9x9"))

        output = self.run_map("m1000")

        self.assertFailedWith(
            output,
            f"{path} was not mapped, nothing of it was stored: The value '9x9' of "
            "well A2 is not a number.",
        )
        self.assertFalse(Measurement.objects.exists())

    def test_an_unexpected_error_shows_its_type_and_is_logged_with_traceback(self):
        with mock.patch(
            "importer.management.commands.map.EchoMapper.run",
            side_effect=KeyError("DMSO"),
        ):
            with self.assertLogs("API", level="ERROR") as logs:
                output = self.run_map("echo")

        self.assertFailedWith(output, "KeyError: 'DMSO'")
        failure = logs.records[-1]
        self.assertEqual("Command map echo failed", failure.getMessage())
        self.assertIsNotNone(failure.exc_info)

    def test_a_folder_that_does_not_exist(self):
        missing = join(self.folder, "missing")

        output = self.run_map("echo", path=missing)

        self.assertFailedWith(output, f"The folder {missing} does not exist.")

    def test_a_file_instead_of_a_folder(self):
        path = self.write("report.csv", "")

        output = self.run_map("echo", path=path)

        self.assertFailedWith(
            output, f"{path} is a file. Please choose the folder that contains it."
        )
