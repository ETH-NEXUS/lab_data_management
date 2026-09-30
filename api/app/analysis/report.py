"""
Renders a statistics report (api/app/statistics/<analysis type>.qmd) with Quarto,
and packs the report and its result files into one zip.
"""

import os
import re
import shutil
import subprocess
import zipfile

import yaml
from django.conf import settings
from django.core.management.base import CommandError

STATISTICS_FOLDER = os.path.join(settings.BASE_DIR, "statistics")

ANALYSIS_TYPES = ["single", "selectivity"]

# The analysis settings a user may choose, with the values the statistics group uses
DEFAULT_SETTINGS = {
    "geom_cor": "median polish",
    "act_cut": "log10(1.5)",
    "fdr_cut": 0.01,
    "select_cut": 0.3,
    "select_cut_yes": "2 * act_cut",
    "select_cut_no": "act_cut",
    "condi_yes": "",
    "condi_no": "",
}
NUMBER_SETTINGS = ["fdr_cut", "select_cut"]

# The reports run these settings as R code (eval(parse(text = ...))), so only
# simple formulas are accepted, e.g. "log10(1.5)" or "2 * act_cut"
FORMULA_SETTINGS = ["act_cut", "select_cut_yes", "select_cut_no"]
FORMULA = re.compile(r"(?:[0-9.\s*/+\-()]|log10|log2|act_cut)+")

# The biggest screen with plate information on prod (75 plates, 28800 wells) takes
# 106 s. The analysis worker runs one report at a time, so a report that hangs is
# stopped after 15 minutes and does not keep the other analyses waiting for long.
RENDER_TIMEOUT_SECONDS = 15 * 60
# Exit codes of GNU timeout: the time was up, the command (quarto) was not found,
# or it had to be killed (128 + 9: by timeout after TERM was not enough, or by the
# system because it used too much memory)
TIMED_OUT = 124
COMMAND_NOT_FOUND = 127
KILLED = 137
# The last lines of the R output that are shown when a report fails
ERROR_LINES = 15

# Quarto colors its output and prints the progress of the report, e.g. "36/75 [plate_qc]"
COLOR_CODE = re.compile(r"\x1b\[[0-9;]*m")
PROGRESS_LINE = re.compile(r"^\s*\d+/\d+(\s|$)")
# Where R stopped, e.g. "Quitting from single.qmd:73-110 [init]"
FAILED_STEP = re.compile(r"Quitting from (\S+) \[([^\]]*)\]")


def check_settings(analysis_type: str, chosen_settings: dict) -> dict:
    """
    The analysis settings with the default for each one that was not chosen.
    Numbers may come as text from the form.

    {"fdr_cut": "0.05"} -> {"geom_cor": "median polish", "act_cut": "log10(1.5)",
                            "fdr_cut": 0.05, ...}
    """
    if analysis_type not in ANALYSIS_TYPES:
        raise CommandError(f"Unknown analysis type: {analysis_type}")
    if not isinstance(chosen_settings, dict):
        raise CommandError(
            'The analysis settings must be names with values, e.g. {"fdr_cut": 0.05}, '
            f"not {chosen_settings!r}."
        )

    checked: dict = dict(DEFAULT_SETTINGS)
    for name, value in chosen_settings.items():
        if name not in DEFAULT_SETTINGS:
            raise CommandError(f"Unknown analysis setting: {name}")
        if value is not None and value != "":
            checked[name] = value

    for name in NUMBER_SETTINGS:
        try:
            checked[name] = float(checked[name])
        except (TypeError, ValueError):
            raise CommandError(f'{name} must be a number, not "{checked[name]}".')
    for name in FORMULA_SETTINGS:
        if not FORMULA.fullmatch(str(checked[name])):
            raise CommandError(
                f'{name} must be a number or a formula like "log10(1.5)" or '
                f'"2 * act_cut", not "{checked[name]}".'
            )
    if analysis_type == "selectivity" and not (
        checked["condi_yes"] and checked["condi_no"]
    ):
        raise CommandError(
            "A selectivity analysis needs both conditions (condi_yes and condi_no)."
        )
    # The report compares the two conditions; R stops with an unclear error for one
    if analysis_type == "selectivity" and checked["condi_yes"] == checked["condi_no"]:
        raise CommandError(
            "A selectivity analysis compares two different conditions, but both are "
            f'"{checked["condi_yes"]}". Choose two different conditions.'
        )
    return checked


def render_report(analysis_type: str, report_params: dict, folder: str) -> str:
    """
    Renders the report into `folder` and returns the path of the html file.
    The html has its figures inside, so it can be opened on its own.

    report_params example (all parameters of the .qmd):
    {"project": "P1", "screen": "Screen 1", "hts_type": "single",
     "path_data": "/x/main_info.csv", "path_output": "/x/output/", ...}
    """
    # Quarto writes the report next to the .qmd, so every run renders its own copy
    report_file = f"{analysis_type}.qmd"
    shutil.copy(os.path.join(STATISTICS_FOLDER, report_file), folder)
    # A file instead of -P options: yaml keeps text as text and numbers as numbers
    with open(os.path.join(folder, "params.yml"), "w") as params_file:
        yaml.safe_dump(report_params, params_file)

    # Quarto runs R in processes of its own. The timeout of subprocess.run would
    # stop Quarto only and leave R running, so GNU timeout is used: it stops the
    # whole group (TERM first, KILL 30 seconds later)
    command = [
        "timeout",
        "--kill-after=30",
        str(RENDER_TIMEOUT_SECONDS),
        "quarto",
        "render",
        report_file,
        "--to",
        "html",
        "--execute-params",
        "params.yml",
        "-M",
        "embed-resources:true",
    ]
    result = subprocess.run(command, cwd=folder, capture_output=True, text=True)

    if result.returncode == TIMED_OUT:
        raise CommandError(
            f"The R report was stopped after {RENDER_TIMEOUT_SECONDS // 60} minutes; "
            "a report usually takes one or two minutes, so it probably hung. Send "
            "this message to the statistics group."
        )
    if result.returncode == KILLED:
        raise CommandError(
            "The R report was killed: either it was still running after "
            f"{RENDER_TIMEOUT_SECONDS // 60} minutes and did not stop, or the system "
            "stopped it because it used too much memory. If it happens again for this "
            "experiment, the celery-analysis container needs more memory."
        )
    if result.returncode == COMMAND_NOT_FOUND:
        raise CommandError(
            "Quarto is not installed in the celery container, so the R report cannot "
            "run. The Docker image has to be built with ENABLE_R=True."
        )
    if result.returncode != 0:
        # R writes its errors to stderr; Quarto itself may write an error to stdout
        raise CommandError(render_error_text(result.stderr + "\n" + result.stdout))
    return os.path.join(folder, f"{analysis_type}.html")


def render_error_text(quarto_output: str) -> str:
    """
    What went wrong in a failed report: the step of the .qmd where R stopped, and
    the R error, without colors and progress lines. Example:
    'The R report stopped with an error in the step "init" (single.qmd:73-110).
     The data passed the checks of LDM, so this is a problem inside the R script ...
     R error:
     Error:
     ! Could not load one or more required packages'
    """
    lines = []
    for line in COLOR_CODE.sub("", quarto_output).splitlines():
        if line.strip() and not PROGRESS_LINE.match(line):
            lines.append(line.rstrip())

    heading = "The R report stopped with an error."
    error_lines = lines[-ERROR_LINES:]
    for index, line in enumerate(lines):
        failed_step = FAILED_STEP.search(line)
        if not failed_step:
            continue
        where, step = failed_step.groups()
        heading = f'The R report stopped with an error in the step "{step}" ({where}).'
        # The R error is printed right before, from the line that starts with "Error"
        error_start = max(0, index - ERROR_LINES)
        for error_index in range(index - 1, error_start - 1, -1):
            if lines[error_index].startswith("Error"):
                error_start = error_index
                break
        error_lines = lines[error_start:index]
        break

    return "\n".join(
        [
            heading,
            "The data passed the checks of LDM, so this is a problem inside the R "
            "script or a case it does not handle. Send this message to the "
            "statistics group.",
            "R error:",
            *error_lines,
        ]
    )


def pack_results(zip_path: str, report_path: str, output_folder: str) -> None:
    """
    One zip with the report, the parameters and every file the report wrote,
    e.g. report.html, params.yml, DAA_results.tsv, plate_stats.tsv, plate_hmap_raw.png.

    The parameters in the zip are without the paths (path_data, path_output, ...):
    they point into the temporary folder of the run and into the container, which
    do not exist for anyone who opens the zip.
    """
    with open(os.path.join(os.path.dirname(report_path), "params.yml")) as params_file:
        report_params = yaml.safe_load(params_file)
    shown_params = {
        name: value
        for name, value in report_params.items()
        if not name.startswith("path_")
    }

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.write(report_path, "report.html")
        archive.writestr("params.yml", yaml.safe_dump(shown_params))
        for file_name in sorted(os.listdir(output_folder)):
            archive.write(os.path.join(output_folder, file_name), file_name)
