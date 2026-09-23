"""
A command whose process was killed (e.g. out of memory) must not stay "running".

The main process of the worker sends task_failure with a WorkerLostError then;
these tests send the signal the same way.
"""

from billiard.exceptions import WorkerLostError
from celery.signals import task_failure
from django.core.cache import cache
from django.test import SimpleTestCase

from importer.command_output import (
    LOST_PROCESS_MESSAGE,
    add_message,
    read_output,
    register_running_command,
    start_command,
)
from importer.tasks import run_management_command


def send_task_failure(sender, exception, form_data):
    """Sends the signal like the worker: the task arguments are (form_data,)."""
    task_failure.send(
        sender=sender,
        task_id="task_1",
        exception=exception,
        args=(form_data,),
        kwargs={},
        traceback=None,
        einfo=None,
    )


class LostProcessTest(SimpleTestCase):
    def setUp(self):
        cache.clear()
        start_command("room_1")
        register_running_command("room_1", "celery@celery")
        add_message("room_1", "info", "Processing file a.csv...")

    def test_a_command_whose_process_was_killed_failed(self):
        send_task_failure(
            run_management_command,
            WorkerLostError("Worker exited prematurely: signal 9 (SIGKILL)."),
            {"command": "map", "room_name": "room_1"},
        )

        output = read_output("room_1", since=0)
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            [
                {"level": "info", "text": "Processing file a.csv..."},
                {"level": "error", "text": LOST_PROCESS_MESSAGE},
                {"level": "error", "text": "Command failed."},
            ],
            output["messages"],
        )
        self.assertEqual({}, cache.get("running_commands"))

    def test_another_error_is_left_to_the_task(self):
        send_task_failure(
            run_management_command,
            ValueError("handled by the task"),
            {"command": "map", "room_name": "room_1"},
        )

        self.assertEqual("running", read_output("room_1", since=0)["status"])

    def test_another_task_is_not_touched(self):
        class OtherTask:
            name = "inventory.tasks.something_else"

        send_task_failure(
            OtherTask(),
            WorkerLostError("Worker exited prematurely: signal 9 (SIGKILL)."),
            {"room_name": "room_1"},
        )

        self.assertEqual("running", read_output("room_1", since=0)["status"])
