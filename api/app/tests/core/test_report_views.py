"""
The report views need a logged in user and only work inside the notebooks folder.
"""

import tempfile
from os.path import join
from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

POST_VIEWS = [
    "generate_pdf_report",
    "download_pdf_report",
    "download_csv_data",
    "list_files",
    "save_plate_info",
]


class ReportViewsTest(TestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()
        notebooks = override_settings(NOTEBOOKS_ROOT=self.folder)
        notebooks.enable()
        self.addCleanup(notebooks.disable)

    def test_without_login_every_report_view_is_refused(self):
        responses = [self.client.post(reverse(name)) for name in POST_VIEWS]
        responses.append(self.client.get(reverse("prefillPlateInfo")))

        self.assertEqual(
            [403] * (len(POST_VIEWS) + 1),
            [response.status_code for response in responses],
        )

    def test_a_report_outside_the_notebooks_folder_is_not_downloaded(self):
        self.client.force_login(User.objects.create_user("tester"))

        response = self.client.post(
            reverse("download_pdf_report"),
            {"path": "/etc/passwd"},
            content_type="application/json",
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("must be inside", response.json()[0])

    def test_a_report_inside_the_notebooks_folder_is_downloaded(self):
        self.client.force_login(User.objects.create_user("tester"))
        path = join(self.folder, "report.pdf")
        with open(path, "wb") as file:
            file.write(b"a report")

        response = self.client.post(
            reverse("download_pdf_report"),
            {"path": path},
            content_type="application/json",
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual(b"a report", response.content)

    def test_a_notebook_outside_the_notebooks_folder_is_not_run(self):
        self.client.force_login(User.objects.create_user("tester"))

        with mock.patch("core.views.reports.subprocess.run") as run:
            response = self.client.post(
                reverse("generate_pdf_report"),
                {"notebook_path": "/etc/passwd", "experiment": "demo"},
                content_type="application/json",
            )

        self.assertEqual(400, response.status_code)
        run.assert_not_called()
