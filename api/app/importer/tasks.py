"""
Background tasks of the importer. They run in the `celery` container, so a long
command is not stopped by the timeout of the web server.
"""

from celery import shared_task
from celery.signals import worker_ready
from django.core import management
from django.core.management.base import CommandError

from helpers.logger import logger
from importer.command_output import (
    error_text,
    fail_interrupted_commands,
    finish_command,
    register_running_command,
)
from importer.helper import message


@shared_task
def run_management_command(form_data: dict) -> None:
    """
    Runs the map or import command that was started on the management page.
    The page reads the messages and the end status through long_polling.

    Accepted data example:
    {"command": "map", "machine": "echo", "path": "/data/run_1",
     "experiment_name": "Screen 1", "room_name": "12_1726563600000"}
    """
    room_name = form_data.get("room_name")
    register_running_command(room_name)
    try:
        if form_data.get("command") == "map":
            kwargs = {
                "path": form_data.get("path"),
                "mapping_file": form_data.get("mapping_file"),
                "debug": False,
                "experiment_name": form_data.get("experiment_name"),
                "room_name": room_name,
                "measurement_name": form_data.get("measurement_name"),
            }
            # An unknown machine is refused by the command itself
            management.call_command("map", form_data.get("machine"), **kwargs)
        elif form_data.get("command") == "import":
            kwargs = {
                "mapping_file": form_data.get("mapping_file"),
                "input_file": form_data.get("input_file"),
                "debug": False,
                "library_name": form_data.get("library_name") or None,
                "template_name": form_data.get("template_name") or None,
                "plate_barcode": form_data.get("plate_barcode") or None,
                "project_name": form_data.get("project_name") or None,
                "is_control_plate": bool(form_data.get("is_control_plate")),
                "room_name": room_name,
            }
            management.call_command("import", form_data.get("what"), **kwargs)
        else:
            raise CommandError(f"Unknown command: {form_data.get('command')}")
    except Exception as error:
        # An error the command did not handle itself, e.g. an unknown experiment
        message(error_text(error), "error", room_name)
        if not isinstance(error, CommandError):
            logger.exception(f"Command {form_data.get('command')} failed")
    finally:
        finish_command(room_name)


@worker_ready.connect
def fail_commands_of_the_last_worker(**kwargs) -> None:
    """A restarted worker does not continue the commands it was running."""
    fail_interrupted_commands()
