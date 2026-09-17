"""
Tests for running a command from the management page and reading its output.
"""

from datetime import datetime
from os.path import join
from unittest import mock

from django.core.cache import cache

from core.models import Experiment, Measurement, Project
from tests.management.base import ManagementPageTestCase

# A shortened C10 reader file; the header line names the values "Lum"
C10_TXT = "\r\n".join(
    [
        "Plate Number\tPlate 1",
        "Date\t{date}",
        "Time\t12:45:28",
        "Results",
        "Well\tLum",
        "A1\t16727",
        "A2\t1.5E+03",
        "",
    ]
)


class RunCommandTest(ManagementPageTestCase):
    def run_map(self, **form_data):
        """Starts a map command and returns the response of the request."""
        data = {
            "command": "map",
            "machine": "echo",
            "path": self.folder,
            "experiment_name": "",
            "measurement_name": "",
        }
        data.update(form_data)
        return self.start_command(**data)

    def test_a_command_without_files_is_completed_with_a_warning(self):
        response = self.run_map()

        self.assertEqual(200, response.status_code)
        output = self.read_output()
        self.assertEqual("completed", output["status"])
        levels_and_texts = [
            (message["level"], message["text"]) for message in output["messages"]
        ]
        self.assertEqual("warning", levels_and_texts[0][0])
        self.assertIn("No files found that match", levels_and_texts[0][1])
        self.assertEqual(("info", "Command completed."), levels_and_texts[-1])

    def test_an_error_outside_the_command_is_shown_once_and_the_command_failed(self):
        # map raises this error before its own error handling
        self.run_map(experiment_name="No such experiment")

        output = self.read_output()
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            [
                {
                    "level": "error",
                    "text": "No experiment with name 'No such experiment' found in the database.",
                },
                {"level": "error", "text": "Command failed."},
            ],
            output["messages"],
        )

    def test_the_command_is_sent_to_celery_and_the_request_returns_at_once(self):
        with mock.patch("management.views.run_management_command.delay") as delay:
            response = self.run_map()

        self.assertEqual(200, response.status_code)
        delay.assert_called_once()
        self.assertEqual("room_1", delay.call_args.args[0]["room_name"])
        # The command has not run yet: no messages, but the page already sees "running"
        self.assertEqual(
            {"messages": [], "next": 0, "status": "running"}, self.read_output()
        )

    def test_a_finished_command_is_no_longer_listed_as_running(self):
        self.run_map()

        self.assertEqual({}, cache.get("running_commands"))

    def test_an_unknown_command_is_an_error(self):
        self.run_map(command="do_something")

        output = self.read_output()
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            {"level": "error", "text": "Unknown command: do_something"},
            output["messages"][0],
        )

    def test_an_unknown_machine_is_an_error(self):
        self.run_map(machine="pipetting_robot")

        output = self.read_output()
        self.assertEqual("failed", output["status"])
        self.assertIn("invalid choice", self.errors(output)[0])

    def test_the_output_can_be_read_from_a_position(self):
        self.run_map()
        all_messages = self.read_output()["messages"]

        output = self.read_output(since=len(all_messages) - 1)

        self.assertEqual([all_messages[-1]], output["messages"])
        self.assertEqual(len(all_messages), output["next"])

    def test_before_a_command_there_is_no_output(self):
        self.assertEqual(
            {"messages": [], "next": 0, "status": None}, self.read_output()
        )

    def write_c10_file(self, date):
        path = join(self.folder, "241014_125455_241008MP-1_1.txt")
        with open(path, "w", newline="") as file:
            file.write(C10_TXT.format(date=date))
        return path

    def test_c10_values_without_measurement_name_are_labeled_from_the_file(self):
        project = Project.objects.create(name="Project")
        Experiment.objects.create(name="Experiment", project=project)
        self.write_c10_file(date="10/14/2024")

        self.run_map(machine="C10-reader", experiment_name="Experiment")

        self.assertEqual("completed", self.read_output()["status"])
        self.assertEqual(
            [
                ("Lum", 16727.0, datetime(2024, 10, 14, 12, 45, 28)),
                ("Lum", 1500.0, datetime(2024, 10, 14, 12, 45, 28)),
            ],
            list(
                Measurement.objects.order_by("well__position").values_list(
                    "label", "value", "measured_at"
                )
            ),
        )

    def test_a_c10_file_with_an_unknown_date_shows_the_error(self):
        project = Project.objects.create(name="Project")
        Experiment.objects.create(name="Experiment", project=project)
        path = self.write_c10_file(date="14.10.2024")

        self.run_map(machine="C10-reader", experiment_name="Experiment")

        output = self.read_output()
        self.assertEqual("failed", output["status"])
        self.assertIn(
            {
                "level": "error",
                "text": f"{path} was not mapped, nothing of it was stored: Cannot read "
                "the measurement date: date '14.10.2024', time '12:45:28'.",
            },
            output["messages"],
        )
        self.assertFalse(Measurement.objects.exists())
