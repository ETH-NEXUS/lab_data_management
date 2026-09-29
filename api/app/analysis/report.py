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

# A report of a whole screen takes a few minutes; this only stops one that hangs
RENDER_TIMEOUT_SECONDS = 60 * 60
# The last lines of the R output that are shown when a report fails
ERROR_LINES = 15

# Quarto colors its output and prints the progress of the report, e.g. "36/75 [plate_qc]"
COLOR_CODE = re.compile(r"\x1b\[[0-9;]*m")
PROGRESS_LINE = re.compile(r"^\s*\d+/\d+(\s|$)")


def check_settings(analysis_type: str, chosen_settings: dict) -> dict:
    """
    The analysis settings with the default for each one that was not chosen.
    Numbers may come as text from the form.

    {"fdr_cut": "0.05"} -> {"geom_cor": "median polish", "act_cut": "log10(1.5)",
                            "fdr_cut": 0.05, ...}
    """
    if analysis_type not in ANALYSIS_TYPES:
        raise CommandError(f"Unknown analysis type: {analysis_type}")

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
    return checked


def render_report(analysis_type: str, report_params: dict, folder: str) -> str:
    """
    Renders the report into `folder` and returns the path of the html file.
    The html has its figures inside, so it can be opened on its own.
    The full Quarto output is kept in render.log.

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

    command = [
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
    log_path = os.path.join(folder, "render.log")
    try:
        result = subprocess.run(
            command,
            cwd=folder,
            capture_output=True,
            text=True,
            timeout=RENDER_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired:
        raise CommandError(
            f"The R report was stopped after {RENDER_TIMEOUT_SECONDS // 60} minutes."
        )

    with open(log_path, "w") as log_file:
        log_file.write(result.stdout + result.stderr)
    if result.returncode != 0:
        raise CommandError(render_error_text(result.stderr, log_path))
    return os.path.join(folder, f"{analysis_type}.html")


def render_error_text(quarto_output: str, log_path: str) -> str:
    """
    The end of the R output of a failed report, without colors and progress lines.
    Example:
    'The R report failed:
     Error:
     ! Could not load one or more required packages
     Quitting from single.qmd:73-110 [init]
     Execution halted
     The full output is in /vol/web/media/analysis/105/.../render.log'
    """
    lines = []
    for line in COLOR_CODE.sub("", quarto_output).splitlines():
        if line.strip() and not PROGRESS_LINE.match(line):
            lines.append(line.rstrip())
    last_lines = lines[-ERROR_LINES:]
    return "\n".join(
        ["The R report failed:", *last_lines, f"The full output is in {log_path}"]
    )


def pack_results(zip_path: str, report_path: str, output_folder: str) -> None:
    """
    One zip with the report, the parameters and every file the report wrote,
    e.g. report.html, params.yml, DAA_results.tsv, plate_stats.tsv, plate_hmap_raw.png.
    """
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.write(report_path, "report.html")
        archive.write(
            os.path.join(os.path.dirname(report_path), "params.yml"), "params.yml"
        )
        for file_name in sorted(os.listdir(output_folder)):
            archive.write(os.path.join(output_folder, file_name), file_name)
