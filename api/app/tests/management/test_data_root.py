"""
The management page only works with files inside the data folder.
"""

import os
import tempfile
from os.path import join
from unittest import mock

from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse

from tests.management.base import ManagementPageTestCase


class DataRootTest(ManagementPageTestCase):
    def setUp(self):
        super().setUp()
        # A file the page must not touch, outside the data folder of this test
        self.outside_file = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        self.outside_file.write(b"a secret")
        self.outside_file.close()
        self.addCleanup(os.remove, self.outside_file.name)

    def post(self, name, data):
        return self.client.post(reverse(name), data, content_type="application/json")

    def test_a_file_outside_the_data_folder_is_not_read(self):
        response = self.post("get_file_content", {"file_path": self.outside_file.name})

        self.assertEqual(400, response.status_code)
        self.assertIn("must be inside", response.json()[0])

    def test_a_file_outside_the_data_folder_is_not_downloaded(self):
        response = self.post("download_file", {"file_path": self.outside_file.name})

        self.assertEqual(400, response.status_code)

    def test_a_file_outside_the_data_folder_is_not_deleted(self):
        response = self.post("delete_file", {"path": self.outside_file.name})

        self.assertEqual(400, response.status_code)
        with open(self.outside_file.name) as file:
            self.assertEqual("a secret", file.read())

    def test_a_path_that_leaves_the_data_folder_is_refused(self):
        response = self.post(
            "get_file_content", {"file_path": join(self.folder, "..", "escape.txt")}
        )

        self.assertEqual(400, response.status_code)

    def test_a_command_with_a_path_outside_the_data_folder_is_not_started(self):
        with mock.patch("management.views.run_management_command.delay") as delay:
            response = self.post(
                "run_command",
                {
                    "form_data": {
                        "command": "import",
                        "what": "sdf",
                        "input_file": "/etc/passwd",
                        "room_name": "room_1",
                    }
                },
            )

        self.assertEqual(400, response.status_code)
        delay.assert_not_called()

    def test_a_file_inside_the_data_folder_is_read(self):
        path = join(self.folder, "notes.txt")
        with open(path, "w") as file:
            file.write("hello")

        response = self.post("get_file_content", {"file_path": path})

        self.assertEqual({"content": "hello"}, response.json())

    def test_an_uploaded_file_keeps_only_its_name(self):
        # A file name that tries to leave the folder it is uploaded into
        uploaded = SimpleUploadedFile("../../escape.txt", b"hello")

        response = self.client.post(
            reverse("upload_file"),
            {"directory_path": self.folder, "file": uploaded},
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual(join(self.folder, "escape.txt"), response.json()["file_path"])
        self.assertEqual(["escape.txt"], os.listdir(self.folder))

    def test_a_request_without_a_path_says_so(self):
        response = self.post("get_file_content", {})

        self.assertEqual(400, response.status_code)
        self.assertEqual(["No path was given."], response.json())

    def test_a_request_without_form_data_says_so(self):
        response = self.post("run_command", {"form_data": None})

        self.assertEqual(400, response.status_code)
        self.assertEqual(["The request has no form_data."], response.json())
