"""
The management page works with the files of the data folder.

Every path comes from the browser, so it is checked before a file is read,
written, deleted or given to a command.
"""

import os

from django.conf import settings
from rest_framework.exceptions import ValidationError


def data_path(path: str) -> str:
    """
    The path, if it is inside the data folder (settings.MANAGEMENT_DATA_ROOT).

    "/data/wagner/report.csv" -> "/data/wagner/report.csv", while
    "/data/../etc/passwd" and "/etc/passwd" are refused with an error that the
    page shows. A symbolic link inside the data folder is allowed: the lab
    shares are mounted and linked in ways this code does not need to know.
    """
    if not path:
        raise ValidationError("No path was given.")

    root = os.path.abspath(settings.MANAGEMENT_DATA_ROOT)
    full_path = os.path.abspath(path)
    if full_path != root and not full_path.startswith(root + os.sep):
        raise ValidationError(f"The path must be inside {root}: {path}")
    return full_path


def check_command_paths(form_data: dict) -> None:
    """Checks the paths of a command of the management page, e.g. its input file."""
    for key in ("path", "input_file", "mapping_file"):
        if form_data.get(key):
            data_path(form_data[key])
