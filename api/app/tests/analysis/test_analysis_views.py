"""
The analysis is started from the experiment page, and its result zips are
listed and downloaded there. Only for logged in users.
"""

import os
import shutil
import tempfile
import time
from datetime import datetime
from unittest import mock

from django.contrib.auth.models import User
from django.core.cache import cache
from django.test import TestCase
from django.urls import reverse

from analysis import views
from core.models import Experiment, Plate, PlateDimension, PlateInfo, Project
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
        (sent,) = delay.call_args.args
        self.assertEqual(form_data, {k: v for k, v in sent.items() if k != "queued_at"})
        # The worker refuses a run that waited too long, so the start time is sent
        self.assertAlmostEqual(time.time(), sent["queued_at"], delta=60)
        self.assertEqual("running", read_output("3_1727600000000", 0)["status"])

    def test_a_room_that_was_started_already_is_refused(self):
        self.login()
        form_data = {
            "experiment_id": 105,
            "label": "Lum_CTG",
            "analysis_type": "single",
            "room_name": "3_1727600000000",
        }
        with mock.patch.object(views.run_analysis, "delay") as delay:
            self.client.post(
                reverse("start_analysis"), form_data, content_type="application/json"
            )
            second = self.client.post(
                reverse("start_analysis"), form_data, content_type="application/json"
            )

        self.assertEqual(400, second.status_code)
        self.assertEqual(1, delay.call_count)

    def test_the_conditions_of_a_measurement_are_listed(self):
        self.login()
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_2x2", rows=2, cols=2)
        for number, condition in enumerate(["irradiated", "not irradiated", "", "KO"]):
            plate = Plate.objects.create(
                barcode=f"SP_{number}", dimension=dimension, experiment=experiment
            )
            PlateInfo.objects.create(
                plate=plate,
                experiment=experiment,
                lib_plate_barcode="LIB_1",
                label="Fluo" if condition == "KO" else "Lum",
                replicate="1",
                measurement_time=datetime(2025, 5, 16),
                cell_type="",
                condition=condition,
            )

        response = self.client.get(
            reverse("list_analysis_conditions"),
            {"experiment_id": experiment.id, "label": "Lum"},
        )

        self.assertEqual(
            {"conditions": ["irradiated", "not irradiated"]}, response.json()
        )

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
