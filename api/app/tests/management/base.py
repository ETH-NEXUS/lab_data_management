"""
What the tests of the management page share.

A test runs a command the way the page does it: it logs in, posts to
`run_command` and reads the output of the command with `long_polling`.
Every test works in its own folder, which is the data folder for that test.
"""

import shutil
import tempfile
from os.path import join

from django.contrib.auth.models import User
from django.core.cache import cache
from django.test import TestCase, override_settings
from django.urls import reverse

ROOM_NAME = "room_1"


class ManagementPageTestCase(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        cache.clear()
        self.client.force_login(User.objects.create_user("tester"))
        self.folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.folder)
        for setting in (
            override_settings(MANAGEMENT_DATA_ROOT=self.folder),
            override_settings(MEDIA_ROOT=join(self.folder, "media")),
        ):
            setting.enable()
            self.addCleanup(setting.disable)

    def write(self, name, text):
        """Writes a file into the data folder of this test and returns its path."""
        path = join(self.folder, name)
        with open(path, "w") as file:
            file.write(text)
        return path

    def start_command(self, **form_data):
        """
        Starts a command as the page does and returns the response. In the tests
        Celery runs the task right away, so the command is done when this returns.
        """
        data = {"room_name": ROOM_NAME}
        data.update(form_data)
        return self.client.post(
            reverse("run_command"),
            {"form_data": data},
            content_type="application/json",
        )

    def read_output(self, since=0):
        """The messages and the status of the command, as the page reads them."""
        url = reverse("long_polling", args=[ROOM_NAME])
        return self.client.get(f"{url}?since={since}").json()

    def errors(self, output):
        """The texts of the error messages of a command."""
        return [
            message["text"]
            for message in output["messages"]
            if message["level"] == "error"
        ]

    def assertFailedWith(self, output, text):
        """The command failed, and `text` is its only error besides "Command failed."."""
        self.assertEqual("failed", output["status"])
        self.assertEqual([text, "Command failed."], self.errors(output))
