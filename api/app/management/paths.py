"""
The management page works with the files of the data folder.

Every path comes from the browser, so it is checked before a file is read,
written, deleted or given to a command.
"""

from django.conf import settings

from helpers.paths import path_inside


def data_path(path: str) -> str:
    """The path, if it is inside the data folder (settings.MANAGEMENT_DATA_ROOT)."""
    return path_inside(path, settings.MANAGEMENT_DATA_ROOT)


def check_command_paths(form_data: dict) -> None:
    """Checks the paths of a command of the management page, e.g. its input file."""
    for key in ("path", "input_file", "mapping_file"):
        if form_data.get(key):
            data_path(form_data[key])
