"""
The analysis task writes the input files of one measurement label, renders the
report (Quarto is not run here) and packs the results into a zip. Every problem
ends as an error message the page can show.
"""

import csv
import os
import shutil
import tempfile
import zipfile
from datetime import datetime
from unittest import mock

from django.core.cache import cache
from django.test import TestCase

from analysis import tasks
from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDimension,
    PlateInfo,
    Project,
    Well,
    WellType,
)
from importer.command_output import read_output, start_command

ROOM = "7_1727600000000"

# What Quarto prints when R stops, with its colors and progress lines
R_ERROR_OUTPUT = (
    "processing file: single.qmd\n"
    "5/75\n"
    "6/75 [init]\n"
    "\x1b[31mError:\n"
    "! Could not load one or more required packages\x1b[39m\n"
    "\x1b[31mQuitting from single.qmd:73-110 [init]\x1b[39m\n"
    "Execution halted\n"
)


def fake_quarto(command, cwd, **kwargs):
    """Writes what the report would write, instead of running Quarto."""
    with open(os.path.join(cwd, "single.html"), "w") as report:
        report.write("<html>report</html>")
    with open(os.path.join(cwd, "output", "DAA_results.tsv"), "w") as results:
        results.write("hits\n")
    return mock.Mock(returncode=0, stdout="Output created: single.html", stderr="")


class AnalysisTaskTest(TestCase):
    def setUp(self):
        cache.clear()
        self.folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.folder)
        patch_folder = mock.patch.object(tasks, "ANALYSIS_FOLDER", self.folder)
        patch_folder.start()
        self.addCleanup(patch_folder.stop)

        project = Project.objects.create(name="P1")
        self.experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_2x2", rows=2, cols=2)
        plate = Plate.objects.create(
            barcode="SP_1", dimension=dimension, experiment=self.experiment
        )
        # Two control wells and two compound wells
        well_types = [
            WellType.objects.create(name=name, description=name)
            for name in ["N", "P", "C"]
        ]
        self.wells = []
        measured_at = datetime(2025, 5, 16, 10, 0)
        for position, well_type in enumerate(well_types + [well_types[2]]):
            well = Well.objects.create(plate=plate, position=position, type=well_type)
            self.wells.append(well)
            for label in ["Lum", "Fluo"]:
                Measurement.objects.create(
                    well=well, label=label, value=position, measured_at=measured_at
                )
        PlateInfo.objects.create(
            plate=plate,
            experiment=self.experiment,
            lib_plate_barcode="LIB_1",
            label="Lum",
            replicate="1",
            measurement_time=measured_at,
            cell_type="SW620",
            condition="irradiated",
        )

    def run_task(self, label="Lum", quarto=fake_quarto):
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": label,
            "analysis_type": "single",
            "settings": {"fdr_cut": "0.05"},
            "room_name": ROOM,
        }
        with mock.patch("analysis.report.subprocess.run", side_effect=quarto):
            tasks.run_analysis.delay(form_data)
        return read_output(ROOM, 0)

    def run_folder(self):
        experiment_folder = os.path.join(self.folder, str(self.experiment.id))
        (run_name,) = os.listdir(experiment_folder)
        return os.path.join(experiment_folder, run_name)

    def error_texts(self, output):
        return [
            message["text"]
            for message in output["messages"]
            if message["level"] == "error"
        ]

    def test_the_input_files_have_only_the_chosen_label(self):
        self.run_task()

        folder = self.run_folder()
        with open(os.path.join(folder, "main_info.csv")) as main_file:
            main_rows = list(csv.DictReader(main_file))
        with open(os.path.join(folder, "experiment_data.csv")) as meta_file:
            meta_rows = list(csv.DictReader(meta_file))
        self.assertEqual(["Lum"] * 4, [row["measurement"] for row in main_rows])
        # The columns SLmisc.R read.LDM() reads
        read_by_the_report = {
            "unique_identifier",
            "plate",
            "plate_row",
            "plate_column",
            "control",
            "value",
        }
        self.assertLessEqual(read_by_the_report, set(main_rows[0]))
        self.assertEqual(1, len(meta_rows))
        self.assertEqual(
            {
                "measurement_label": "Lum",
                "plate": "SP_1",
                "lib_plate_barcode": "LIB_1",
                "cell_type": "SW620",
                "condition": "irradiated",
            },
            {
                name: meta_rows[0][name]
                for name in [
                    "measurement_label",
                    "plate",
                    "lib_plate_barcode",
                    "cell_type",
                    "condition",
                ]
            },
        )

    def test_the_zip_has_the_report_the_parameters_and_the_results(self):
        output = self.run_task()

        folder = self.run_folder()
        zip_name = os.path.basename(folder) + ".zip"
        with zipfile.ZipFile(os.path.join(folder, zip_name)) as archive:
            names = archive.namelist()
            params = archive.read("params.yml").decode()
        self.assertEqual(["report.html", "params.yml", "DAA_results.tsv"], names)
        self.assertIn("fdr_cut: 0.05", params)
        self.assertEqual("completed", output["status"])
        self.assertEqual("success", output["messages"][-2]["level"])

    def test_a_label_without_plate_information_names_the_labels_that_have_it(self):
        output = self.run_task(label="Fluo")

        self.assertEqual("failed", output["status"])
        (error,) = self.error_texts(output)[:1]
        self.assertIn('has no plate information for the measurement "Fluo"', error)
        self.assertIn('with "add experiment data"', error)
        self.assertIn("Plate information exists for: Lum.", error)

    def test_without_positive_controls_the_report_does_not_start(self):
        WellType.objects.filter(name="P").update(name="X")

        with mock.patch("analysis.report.subprocess.run") as quarto:
            output = self.run_task(quarto=quarto)

        quarto.assert_not_called()
        self.assertEqual("failed", output["status"])
        self.assertIn(
            "No plate of this measurement has positive (P) control wells",
            self.error_texts(output)[0],
        )

    def test_a_selectivity_with_an_unknown_condition_names_the_known_ones(self):
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": "Lum",
            "analysis_type": "selectivity",
            "settings": {"condi_yes": "irradiated", "condi_no": "not irradiated"},
            "room_name": ROOM,
        }
        with mock.patch("analysis.report.subprocess.run") as quarto:
            tasks.run_analysis.delay(form_data)

        quarto.assert_not_called()
        self.assertIn(
            "Not in the plate information of this measurement: 'not irradiated'. "
            "Its conditions are: 'irradiated'.",
            self.error_texts(read_output(ROOM, 0))[0],
        )

    def test_an_r_error_names_the_failed_step_and_the_log(self):
        def failing_quarto(command, cwd, **kwargs):
            return mock.Mock(returncode=1, stdout="", stderr=R_ERROR_OUTPUT)

        output = self.run_task(quarto=failing_quarto)

        log_path = os.path.join(self.run_folder(), "render.log")
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            'The R report failed in the step "init" (single.qmd:73-110):\n'
            "Error:\n"
            "! Could not load one or more required packages\n"
            f"The full output is in {log_path}",
            self.error_texts(output)[0],
        )
        self.assertTrue(os.path.exists(log_path))

    def test_an_unexpected_error_says_where_to_find_the_details(self):
        with mock.patch.object(tasks, "pack_results", side_effect=KeyError("x")):
            output = self.run_task()

        self.assertEqual("failed", output["status"])
        self.assertEqual(
            "The analysis stopped because of an unexpected error in LDM: "
            "KeyError: 'x'. The details are in the log of the celery container.",
            self.error_texts(output)[0],
        )
