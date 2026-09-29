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

    def run_task(self, label="Lum", quarto=fake_quarto, **controls):
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": label,
            "analysis_type": "single",
            "settings": {"fdr_cut": "0.05"},
            "room_name": ROOM,
            **controls,
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
        self.assertEqual(
            'The experiment "Screen 1" has no plate information for the measurement '
            '"Fluo". The report needs it for every plate (library plate, replicate, '
            "cell type, condition). Add it on the experiment page with "
            '"add experiment data" and save it. Plate information exists for: "Lum".',
            self.error_texts(output)[0],
        )

    def test_a_label_without_measurements_names_the_labels_there_are(self):
        PlateInfo.objects.update(label="Other")

        output = self.run_task(label="Other")

        self.assertEqual(
            'The experiment "Screen 1" has no measurements "Other". '
            'Its measurements are: "Fluo", "Lum".',
            self.error_texts(output)[0],
        )

    def test_a_plate_without_controls_gets_a_warning(self):
        # A second plate with compound wells only; SP_1 has its controls
        plate = Plate.objects.create(
            barcode="SP_2",
            dimension=self.wells[0].plate.dimension,
            experiment=self.experiment,
        )
        well = Well.objects.create(plate=plate, position=0, type=self.wells[2].type)
        Measurement.objects.create(
            well=well, label="Lum", value=1, measured_at=datetime(2025, 5, 16, 10, 0)
        )

        output = self.run_task()

        warnings = [m["text"] for m in output["messages"] if m["level"] == "warning"]
        self.assertEqual(
            [
                "These plates have measurements but no plate information, so the "
                "report has no library plate, replicate, cell type or condition for "
                'them: SP_2. Add it with "add experiment data".',
                'Plate SP_2 has no positive control wells ("P") and no negative '
                'control wells ("N"), so the report cannot normalize it and leaves '
                "it out of the results.",
            ],
            warnings,
        )
        self.assertEqual("completed", output["status"])

    def test_numbered_controls_are_written_as_p_and_n(self):
        WellType.objects.filter(name="P").update(name="P1")
        WellType.objects.filter(name="N").update(name="N1")

        output = self.run_task(positive_control="P1", negative_control="N1")

        with open(os.path.join(self.run_folder(), "main_info.csv")) as main_file:
            controls = [row["control"] for row in csv.DictReader(main_file)]
        self.assertEqual(["N", "P", "C", "C"], controls)
        self.assertEqual("completed", output["status"])

    def test_a_p_that_was_not_chosen_is_not_a_control_for_r(self):
        WellType.objects.filter(name="C").update(name="P1")

        self.run_task(positive_control="P1", negative_control="N")

        with open(os.path.join(self.run_folder(), "main_info.csv")) as main_file:
            controls = [row["control"] for row in csv.DictReader(main_file)]
        self.assertEqual(["N", "P (not chosen)", "P", "P"], controls)

    def test_a_chosen_control_that_is_not_there_names_the_well_types(self):
        with mock.patch("analysis.report.subprocess.run") as quarto:
            output = self.run_task(quarto=quarto, positive_control="P1")

        quarto.assert_not_called()
        self.assertEqual(
            'The positive control "P1" is not a well type of this measurement. '
            'Its well types are: "C", "N", "P". Choose the control wells in the '
            "analysis window.",
            self.error_texts(output)[0],
        )

    def test_without_both_controls_on_any_plate_r_is_not_started(self):
        # SP_1 keeps its P well only, SP_2 gets the N well
        Well.objects.filter(pk=self.wells[0].pk).update(type=self.wells[2].type)
        plate = Plate.objects.create(
            barcode="SP_2",
            dimension=self.wells[0].plate.dimension,
            experiment=self.experiment,
        )
        well = Well.objects.create(plate=plate, position=0, type=self.wells[0].type)
        Measurement.objects.create(
            well=well, label="Lum", value=1, measured_at=datetime(2025, 5, 16, 10, 0)
        )

        with mock.patch("analysis.report.subprocess.run") as quarto:
            output = self.run_task(quarto=quarto)

        quarto.assert_not_called()
        self.assertEqual(
            'No plate of this measurement has both the positive ("P") and the '
            'negative ("N") control wells, so the report cannot normalize any plate. '
            "Check the controls chosen in the analysis window and the well types of "
            "the plate layout.\n"
            'Plate SP_1 has no negative control wells ("N"), so the report cannot '
            "normalize it and leaves it out of the results.\n"
            'Plate SP_2 has no positive control wells ("P"), so the report cannot '
            "normalize it and leaves it out of the results.",
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
        self.assertEqual(
            "A selectivity analysis compares two conditions of the plate information, "
            'but these are not in it: "not irradiated". The conditions of this '
            'measurement are: "irradiated". Type two of them exactly like this, or '
            'correct the conditions with "add experiment data".',
            self.error_texts(read_output(ROOM, 0))[0],
        )

    def test_an_r_error_names_the_step_the_r_error_and_the_log(self):
        def failing_quarto(command, cwd, **kwargs):
            return mock.Mock(returncode=1, stdout="", stderr=R_ERROR_OUTPUT)

        output = self.run_task(quarto=failing_quarto)

        log_path = os.path.join(self.run_folder(), "render.log")
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            'The R report stopped with an error in the step "init" '
            "(single.qmd:73-110).\n"
            "The data passed the checks of LDM, so this is a problem inside the R "
            "script or a case it does not handle. Send the full output to the "
            "statistics group.\n"
            "R error:\n"
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
