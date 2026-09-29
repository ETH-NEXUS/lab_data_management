"""
The analysis is started from the experiment page, and its result zips are
listed and downloaded there. Only for logged in users.
"""

import os
import shutil
import tempfile
from unittest import mock

from django.contrib.auth.models import User
from django.core.cache import cache
from django.test import TestCase
from django.urls import reverse

from analysis import views
from importer.command_output import read_output

RUN_NAME = "20260929-125825_single_Lum_CTG"


class AnalysisViewsTest(TestCase):
    def setUp(self):
        cache.clear()
        self.folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.folder)
        patch_folder = mock.patch.object(views, "ANALYSIS_FOLDER", self.folder)
        patch_folder.start()
        self.addCleanup(patch_folder.stop)
        # The zip of one finished run, and a file that is not a result
        os.makedirs(os.path.join(self.folder, "105"))
        with open(os.path.join(self.folder, "105", f"{RUN_NAME}.zip"), "wb") as f:
            f.write(b"a zip")
        with open(os.path.join(self.folder, "105", "notes.txt"), "w") as f:
            f.write("not a result")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def test_without_login_nothing_is_started_or_shown(self):
        responses = [
            self.client.post(reverse("start_analysis")),
            self.client.get(reverse("list_analysis_results"), {"experiment_id": 105}),
            self.client.get(reverse("download_analysis_result")),
        ]

        self.assertEqual([403, 403, 403], [r.status_code for r in responses])

    def test_start_runs_the_task_and_the_output_is_running(self):
        self.login()
        form_data = {
            "experiment_id": 105,
            "label": "Lum_CTG",
            "analysis_type": "single",
            "room_name": "3_1727600000000",
        }

        with mock.patch.object(views.run_analysis, "delay") as delay:
            response = self.client.post(
                reverse("start_analysis"), form_data, content_type="application/json"
            )

        self.assertEqual(200, response.status_code)
        delay.assert_called_once_with(form_data)
        self.assertEqual("running", read_output("3_1727600000000", 0)["status"])

    def test_only_finished_runs_are_listed(self):
        self.login()

        response = self.client.get(
            reverse("list_analysis_results"), {"experiment_id": 105}
        )

        self.assertEqual({"results": [f"{RUN_NAME}.zip"]}, response.json())

    def test_a_result_is_downloaded(self):
        self.login()

        response = self.client.get(
            reverse("download_analysis_result"),
            {"experiment_id": 105, "name": f"{RUN_NAME}.zip"},
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual(b"a zip", b"".join(response.streaming_content))

    def test_a_path_is_not_a_result_name(self):
        self.login()

        response = self.client.get(
            reverse("download_analysis_result"),
            {"experiment_id": 105, "name": "../../../etc/passwd.zip"},
        )

        self.assertEqual(400, response.status_code)
