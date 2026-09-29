"""
The statistical analysis of an experiment. It runs in the `celery` container,
because a report of a whole screen takes several minutes.
"""

import os
import re
import shutil
import tempfile
import time
from datetime import datetime

from celery import shared_task
from celery.exceptions import WorkerLostError
from celery.signals import task_failure
from django.conf import settings
from django.core.management.base import CommandError

from analysis.input_files import write_input_files
from analysis.report import (
    STATISTICS_FOLDER,
    check_settings,
    pack_results,
    render_report,
)
from core.models import Experiment
from helpers.logger import logger
from importer.command_output import (
    error_text,
    fail_lost_command,
    finish_command,
    register_running_command,
)
from importer.helper import message

# The zip of every finished run: <MEDIA_ROOT>/analysis/<experiment id>/<run name>.zip
ANALYSIS_FOLDER = os.path.join(settings.MEDIA_ROOT, "analysis")

# A run starts within seconds, or after one or two other analyses (about 2 minutes
# each). One that waited longer, is not started: the worker was not running. The
# page stops waiting a minute later (START_GIVE_UP_MS in ui/app/stores/analysis.ts).
MAX_WAITING_SECONDS = 10 * 60

LOST_PROCESS_MESSAGE = (
    "The analysis was stopped, because the process that ran it was killed "
    "(for example, the R report used too much memory). No result was saved; "
    "you can start it again. If it happens again, the experiment may be too big "
    "for the memory of the celery container."
)


@shared_task(bind=True)
def run_analysis(self, form_data: dict) -> None:
    """
    Makes the report of one experiment and one measurement label. The messages
    and the end status are read through long_polling, as for the import commands.

    Accepted data example:
    {"experiment_id": 105, "label": "Lum_CTG", "analysis_type": "single",
     "settings": {"act_cut": "log10(1.5)", "fdr_cut": 0.01},
     "positive_control": "P1", "negative_control": "N1",
     "room_name": "105_1727600000000"}
    """
    room_name = form_data.get("room_name")
    # A restarted worker then ends this run as failed, like an import command
    register_running_command(room_name, self.request.hostname or "unknown worker")
    try:
        check_waiting_time(form_data)
        message("Running the statistical analysis", "info", room_name)
        zip_path = make_analysis(form_data, room_name)
        message(
            f"The analysis is done: {os.path.basename(zip_path)}", "success", room_name
        )
    except CommandError as error:
        # Written for the user: what is wrong and what to do
        message(str(error), "error", room_name)
    except Exception as error:
        message(
            f"The analysis stopped because of an unexpected error in LDM: "
            f"{error_text(error)}. The details are in the log of the celery container.",
            "error",
            room_name,
        )
        logger.exception(f"Analysis failed: {form_data}")
    finally:
        finish_command(room_name)


@task_failure.connect
def fail_the_analysis_of_a_lost_process(
    sender=None, exception=None, args=None, **kwargs
) -> None:
    """
    The process that ran an analysis was killed, e.g. because R used too much
    memory. The `finally` of the task then never runs, and the page would wait
    forever. The main process of the worker still gets this signal and ends the
    analysis instead (the same as fail_the_command_of_a_lost_process of the importer).
    """
    if sender is None or sender.name != run_analysis.name:
        return
    if not isinstance(exception, WorkerLostError):
        return
    # The task was started with run_analysis.delay(form_data)
    form_data = args[0] if args else {}
    fail_lost_command(form_data.get("room_name"), LOST_PROCESS_MESSAGE)


def check_waiting_time(form_data: dict) -> None:
    """
    An analysis that waited too long in the queue (the analysis worker was not
    running) is not started any more: the page stops waiting for it a minute
    later, and nobody expects its result hours later.
    `queued_at` is set by the start view, e.g. 1727600000.5 (seconds).
    """
    waited_seconds = time.time() - float(form_data.get("queued_at") or time.time())
    if waited_seconds > MAX_WAITING_SECONDS:
        raise CommandError(
            f"The analysis was not started, because it waited {waited_seconds // 60:.0f} "
            "minutes for the analysis worker (container celery-analysis), which was "
            "probably not running. Start it again."
        )


def make_analysis(form_data: dict, room_name: str | None) -> str:
    """Writes the input files, renders the report and returns the path of the zip."""
    analysis_type = form_data.get("analysis_type") or ""
    label = form_data.get("label") or ""
    chosen_settings = check_settings(analysis_type, form_data.get("settings") or {})
    try:
        experiment = Experiment.objects.get(pk=form_data.get("experiment_id"))
    except Experiment.DoesNotExist:
        raise CommandError(f"Experiment {form_data.get('experiment_id')} not found.")

    # e.g. "20260929-101500_single_Lum_CTG"; the label may have spaces or slashes
    safe_label = re.sub(r"[^A-Za-z0-9_-]+", "_", label)
    run_name = f"{datetime.now():%Y%m%d-%H%M%S}_{analysis_type}_{safe_label}"

    # The well types of the controls; the reports know them as "P" and "N"
    controls = {
        "positive": form_data.get("positive_control") or "P",
        "negative": form_data.get("negative_control") or "N",
    }
    conditions = []
    if analysis_type == "selectivity":
        conditions = [chosen_settings["condi_yes"], chosen_settings["condi_no"]]

    # The run works in a temporary folder, which is deleted at the end, also when
    # the run fails. Only the finished zip is kept in the media folder.
    with tempfile.TemporaryDirectory() as folder:
        output_folder = os.path.join(folder, "output")
        os.makedirs(output_folder)

        message(
            f'Step 1 of 3: collecting the data of "{experiment.name}", measurement '
            f'"{label}", positive control "{controls["positive"]}", negative control '
            f'"{controls["negative"]}"',
            "info",
            room_name,
        )
        input_paths, warnings = write_input_files(
            experiment, label, folder, conditions, controls
        )
        for warning in warnings:
            message(warning, "warning", room_name)

        message(
            "Step 2 of 3: making the R report (this takes a few minutes)",
            "info",
            room_name,
        )
        report_params = {
            "project": experiment.project.name,
            "screen": experiment.name,
            "hts_type": analysis_type,
            **chosen_settings,
            **input_paths,
            # The reports add file names to these two paths, so they end with "/"
            "path_output": output_folder + "/",
            "path_SLmisc": STATISTICS_FOLDER + "/",
        }
        report_path = render_report(analysis_type, report_params, folder)

        message(
            "Step 3 of 3: packing the report and the result files", "info", room_name
        )
        temporary_zip_path = os.path.join(folder, f"{run_name}.zip")
        pack_results(temporary_zip_path, report_path, output_folder)

        # /tmp and the media folder are different disks, so the zip is copied
        # under a name the page does not list (".part") and then renamed: a rename
        # in one folder happens at once, so the page never lists a half copied zip
        experiment_folder = os.path.join(ANALYSIS_FOLDER, str(experiment.id))
        os.makedirs(experiment_folder, exist_ok=True)
        zip_path = unused_zip_path(experiment_folder, run_name)
        partial_path = zip_path + ".part"
        shutil.copyfile(temporary_zip_path, partial_path)
        os.replace(partial_path, zip_path)
    return zip_path


def unused_zip_path(folder: str, run_name: str) -> str:
    """
    The path of the zip of a run, with a number added when two runs of the same
    measurement were started in the same second, so none replaces the other.

    "/x/20260929-101500_single_Lum.zip", or "/x/20260929-101500_single_Lum_2.zip"
    when the first one exists
    """
    zip_path = os.path.join(folder, f"{run_name}.zip")
    number = 2
    while os.path.exists(zip_path):
        zip_path = os.path.join(folder, f"{run_name}_{number}.zip")
        number += 1
    return zip_path
