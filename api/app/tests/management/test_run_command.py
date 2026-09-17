"""
Tests for running a command from the management page and reading its output.
"""

import shutil
import tempfile

from django.core.cache import cache
from django.test import TestCase
from django.urls import reverse


class RunCommandTest(TestCase):
    fixtures = ["plate_dimensions"]

    def setUp(self):
        cache.clear()
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def run_map(self, **form_data):
        data = {
            "command": "map",
            "machine": "echo",
            "path": self.folder,
            "experiment_name": "",
            "measurement_name": "",
            "room_name": "room_1",
        }
        data.update(form_data)
        return self.client.post(
            reverse("run_command"),
            {"form_data": data},
            content_type="application/json",
        )

    def read_output(self, since=0):
        url = reverse("long_polling", args=["room_1"])
        return self.client.get(f"{url}?since={since}").json()

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
