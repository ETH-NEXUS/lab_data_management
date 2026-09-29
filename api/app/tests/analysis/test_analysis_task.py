"""
The analysis task writes the input files of one measurement label, renders the
report (Quarto is not run here) and packs the results into a zip. Every problem
ends as an error message the page can show.
"""

import csv
import os
import shutil
import tempfile
import time
import zipfile
from datetime import datetime
from unittest import mock

import pandas as pd
from celery.exceptions import WorkerLostError
from django.core.cache import cache
from django.db import transaction
from django.test import TestCase

from analysis import tasks
from analysis.input_files import MAIN_COLUMNS
from compoundlib.models import Compound
from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDimension,
    PlateInfo,
    Project,
    Well,
    WellCompound,
    WellType,
)
from importer.command_output import read_output, start_command
from ldm.ldm import get_experiment_measurements

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


class AnalysisTaskTest(TestCase):
    def setUp(self):
        cache.clear()
        self.folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.folder)
        # Copies of the input files: the folder of a run is deleted at its end
        self.inputs_folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.inputs_folder)
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

    def fake_quarto(self, command, cwd, **kwargs):
        """
        Writes what the report would write, instead of running Quarto, and keeps
        a copy of the input files for the test.
        """
        for name in ["main_info.csv", "experiment_data.csv"]:
            shutil.copy(os.path.join(cwd, name), self.inputs_folder)
        with open(os.path.join(cwd, "single.html"), "w") as report:
            report.write("<html>report</html>")
        with open(os.path.join(cwd, "output", "DAA_results.tsv"), "w") as results:
            results.write("hits\n")
        return mock.Mock(returncode=0, stdout="Output created: single.html", stderr="")

    def run_task(self, label="Lum", quarto=None, **controls):
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": label,
            "analysis_type": "single",
            "settings": {"fdr_cut": "0.05"},
            "room_name": ROOM,
            **controls,
        }
        quarto = quarto or self.fake_quarto
        with mock.patch("analysis.report.subprocess.run", side_effect=quarto):
            tasks.run_analysis.delay(form_data)
        return read_output(ROOM, 0)

    def read_input(self, name):
        """The rows of one input file of the last run, e.g. read_input("main_info.csv")."""
        with open(os.path.join(self.inputs_folder, name)) as input_file:
            return list(csv.DictReader(input_file))

    def saved_files(self):
        """The files the runs of this experiment left in the media folder."""
        experiment_folder = os.path.join(self.folder, str(self.experiment.id))
        if not os.path.isdir(experiment_folder):
            return []
        return os.listdir(experiment_folder)

    def error_texts(self, output):
        return [
            message["text"]
            for message in output["messages"]
            if message["level"] == "error"
        ]

    def test_the_input_files_have_only_the_chosen_label(self):
        self.run_task()

        main_rows = self.read_input("main_info.csv")
        meta_rows = self.read_input("experiment_data.csv")
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

    def test_only_the_zip_is_kept_with_the_report_the_parameters_and_the_results(self):
        output = self.run_task()

        (zip_name,) = self.saved_files()
        zip_path = os.path.join(self.folder, str(self.experiment.id), zip_name)
        with zipfile.ZipFile(zip_path) as archive:
            names = archive.namelist()
            params = archive.read("params.yml").decode()
        self.assertTrue(zip_name.endswith("_single_Lum.zip"))
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

        controls = [row["control"] for row in self.read_input("main_info.csv")]
        self.assertEqual(["N", "P", "C", "C"], controls)
        self.assertEqual("completed", output["status"])

    def test_a_p_that_was_not_chosen_is_not_a_control_for_r(self):
        WellType.objects.filter(name="C").update(name="P1")

        self.run_task(positive_control="P1", negative_control="N")

        controls = [row["control"] for row in self.read_input("main_info.csv")]
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

    def test_a_selectivity_without_any_condition_says_to_fill_them_in(self):
        PlateInfo.objects.update(condition="")
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": "Lum",
            "analysis_type": "selectivity",
            "settings": {"condi_yes": "irradiated", "condi_no": "not irradiated"},
            "room_name": ROOM,
        }
        with mock.patch("analysis.report.subprocess.run"):
            tasks.run_analysis.delay(form_data)

        self.assertEqual(
            "A selectivity analysis compares two conditions of the plate information, "
            "but no plate of this measurement has a condition yet. Fill in the column "
            '"Condition" with "add experiment data" and save it.',
            self.error_texts(read_output(ROOM, 0))[0],
        )

    def test_an_r_error_names_the_step_and_the_r_error_and_leaves_no_files(self):
        def failing_quarto(command, cwd, **kwargs):
            return mock.Mock(returncode=1, stdout="", stderr=R_ERROR_OUTPUT)

        output = self.run_task(quarto=failing_quarto)

        self.assertEqual("failed", output["status"])
        self.assertEqual(
            'The R report stopped with an error in the step "init" '
            "(single.qmd:73-110).\n"
            "The data passed the checks of LDM, so this is a problem inside the R "
            "script or a case it does not handle. Send this message to the "
            "statistics group.\n"
            "R error:\n"
            "Error:\n"
            "! Could not load one or more required packages",
            self.error_texts(output)[0],
        )
        self.assertEqual([], self.saved_files())

    def test_main_info_from_the_chemical_export_is_the_main_export(self):
        # A compound with library data in one well, as in a real screen
        compound = Compound.objects.create(name="Cpd 1", data={"ID": "X1", "MW": 300.5})
        WellCompound.objects.create(well=self.wells[2], compound=compound, amount=10)

        main = get_experiment_measurements(
            "Screen 1", "Lum", "main", csv=True, experiment_id=self.experiment.id
        )
        chemical = get_experiment_measurements(
            "Screen 1", "Lum", "chemical", csv=True, experiment_id=self.experiment.id
        )

        pd.testing.assert_frame_equal(main, chemical[MAIN_COLUMNS])

    def test_a_file_mapped_twice_is_used_once_with_a_warning(self):
        # The same values once more, one day later
        for well in self.wells:
            Measurement.objects.create(
                well=well,
                label="Lum",
                value=well.position,
                measured_at=datetime(2025, 5, 17, 10, 0),
            )

        output = self.run_task()

        self.assertEqual(4, len(self.read_input("main_info.csv")))
        warnings = [m["text"] for m in output["messages"] if m["level"] == "warning"]
        self.assertEqual(
            [
                'The measurement "Lum" was mapped more than once: 4 readings are copies '
                "of another reading of the same well with the same value. Each well is "
                "used once in the report."
            ],
            warnings,
        )
        self.assertEqual("completed", output["status"])

    def test_two_readings_with_different_values_stop_the_analysis(self):
        Measurement.objects.create(
            well=self.wells[0],
            label="Lum",
            value=99,
            measured_at=datetime(2025, 5, 17, 10, 0),
        )

        with mock.patch("analysis.report.subprocess.run") as quarto:
            output = self.run_task(quarto=quarto)

        quarto.assert_not_called()
        self.assertEqual(
            'The measurement "Lum" has more than one reading with different values on '
            "1 well(s) (e.g. well A1 of plate SP_1). The report can use one reading per "
            "well only: map the readings with different measurement names, or keep "
            "only the reading that belongs to the analysis.",
            self.error_texts(output)[0],
        )

    def test_the_results_are_deleted_with_the_experiment(self):
        self.run_task()
        self.assertEqual(1, len(self.saved_files()))
        experiment_folder = os.path.join(self.folder, str(self.experiment.id))

        with self.captureOnCommitCallbacks(execute=True):
            self.experiment.delete()

        self.assertFalse(os.path.exists(experiment_folder))

    def test_the_results_stay_when_the_deletion_is_rolled_back(self):
        self.run_task()
        experiment_folder = os.path.join(self.folder, str(self.experiment.id))

        with self.assertRaises(RuntimeError):
            with transaction.atomic():
                self.experiment.delete()
                raise RuntimeError("a later step of the same request fails")

        self.assertTrue(os.path.exists(experiment_folder))
        self.assertTrue(Experiment.objects.filter(name="Screen 1").exists())

    def test_a_killed_report_says_why_it_may_have_been_killed(self):
        def killed_quarto(command, cwd, **kwargs):
            # GNU timeout ends with 128 + 9 when the command was killed
            return mock.Mock(returncode=137, stdout="", stderr="")

        output = self.run_task(quarto=killed_quarto)

        self.assertIn("The R report was killed", self.error_texts(output)[0])
        self.assertIn("too much memory", self.error_texts(output)[0])

    def test_a_run_that_waited_too_long_in_the_queue_is_not_started(self):
        start_command(ROOM)
        form_data = {
            "experiment_id": self.experiment.id,
            "label": "Lum",
            "analysis_type": "single",
            "room_name": ROOM,
            "queued_at": time.time() - 11 * 60,
        }
        with mock.patch("analysis.report.subprocess.run") as quarto:
            tasks.run_analysis.delay(form_data)

        quarto.assert_not_called()
        output = read_output(ROOM, 0)
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            "The analysis was not started, because it waited 11 minutes for the "
            "analysis worker (container celery-analysis), which was probably not "
            "running. Start it again.",
            self.error_texts(output)[0],
        )
        self.assertEqual([], self.saved_files())

    def test_a_plate_whose_first_well_was_not_read_is_exported(self):
        # The first well (an N control) was not read; another well is the N control
        Measurement.objects.filter(well=self.wells[0]).delete()
        Well.objects.filter(pk=self.wells[3].pk).update(type=self.wells[0].type)

        self.run_task()

        plates = {row["plate"] for row in self.read_input("main_info.csv")}
        self.assertEqual({"SP_1"}, plates)
        self.assertEqual(3, len(self.read_input("main_info.csv")))

    def test_a_report_that_takes_too_long_says_so(self):
        def stopped_quarto(command, cwd, **kwargs):
            # GNU timeout ends with 124 when the time was up
            return mock.Mock(returncode=124, stdout="", stderr="")

        output = self.run_task(quarto=stopped_quarto)

        self.assertEqual(
            "The R report was stopped after 15 minutes; a report usually takes one or "
            "two minutes, so it probably hung. Send this message to the statistics "
            "group.",
            self.error_texts(output)[0],
        )

    def test_without_quarto_the_error_says_how_to_get_it(self):
        def missing_quarto(command, cwd, **kwargs):
            # GNU timeout ends with 127 when the command does not exist
            return mock.Mock(returncode=127, stdout="", stderr="")

        output = self.run_task(quarto=missing_quarto)

        self.assertIn("Quarto is not installed", self.error_texts(output)[0])

    def test_a_second_zip_of_the_same_second_gets_a_number(self):
        folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, folder)
        first = tasks.unused_zip_path(folder, "20260929-101500_single_Lum")
        open(first, "w").close()

        second = tasks.unused_zip_path(folder, "20260929-101500_single_Lum")

        self.assertEqual(
            ["20260929-101500_single_Lum.zip", "20260929-101500_single_Lum_2.zip"],
            [os.path.basename(first), os.path.basename(second)],
        )

    def test_an_experiment_with_the_same_name_in_another_project_is_not_mixed_in(self):
        other_project = Project.objects.create(name="P2")
        other_experiment = Experiment.objects.create(
            name="Screen 1", project=other_project
        )
        plate = Plate.objects.create(
            barcode="OTHER_1",
            dimension=self.wells[0].plate.dimension,
            experiment=other_experiment,
        )
        well = Well.objects.create(plate=plate, position=0, type=self.wells[0].type)
        Measurement.objects.create(
            well=well, label="Lum", value=1, measured_at=datetime(2025, 5, 16, 10, 0)
        )

        self.run_task()

        plates = {row["plate"] for row in self.read_input("main_info.csv")}
        self.assertEqual({"SP_1"}, plates)

    def test_a_killed_process_ends_the_analysis_as_failed(self):
        start_command(ROOM)

        tasks.fail_the_analysis_of_a_lost_process(
            sender=tasks.run_analysis,
            exception=WorkerLostError(),
            args=[{"room_name": ROOM}],
        )

        output = read_output(ROOM, 0)
        self.assertEqual("failed", output["status"])
        self.assertEqual(tasks.LOST_PROCESS_MESSAGE, self.error_texts(output)[0])

    def test_another_error_of_the_task_is_left_to_the_task(self):
        start_command(ROOM)

        tasks.fail_the_analysis_of_a_lost_process(
            sender=tasks.run_analysis,
            exception=KeyError("x"),
            args=[{"room_name": ROOM}],
        )

        self.assertEqual("running", read_output(ROOM, 0)["status"])

    def test_an_unexpected_error_says_where_to_find_the_details(self):
        with mock.patch.object(tasks, "pack_results", side_effect=KeyError("x")):
            output = self.run_task()

        self.assertEqual("failed", output["status"])
        self.assertEqual(
            "The analysis stopped because of an unexpected error in LDM: "
            "KeyError: 'x'. The details are in the log of the celery container.",
            self.error_texts(output)[0],
        )
