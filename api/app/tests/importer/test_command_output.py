"""
Tests for the output of a command started from the management page.
"""

from django.core.cache import cache
from django.core.management.base import CommandError
from django.test import SimpleTestCase

from importer.command_output import (
    INTERRUPTED_MESSAGE,
    MAX_MESSAGES,
    TOO_MANY_MESSAGES,
    add_message,
    error_text,
    fail_interrupted_commands,
    finish_command,
    read_output,
    register_running_command,
    start_command,
)


class CommandOutputTest(SimpleTestCase):
    def setUp(self):
        cache.clear()

    def test_messages_are_kept_with_their_level(self):
        start_command("room_1")
        add_message("room_1", "info", "Processing file a.csv...")
        add_message("room_1", "warning", "No files found.")

        self.assertEqual(
            {
                "messages": [
                    {"level": "info", "text": "Processing file a.csv..."},
                    {"level": "warning", "text": "No files found."},
                ],
                "next": 2,
                "status": "running",
            },
            read_output("room_1", since=0),
        )

    def test_only_new_messages_are_read(self):
        start_command("room_1")
        add_message("room_1", "info", "first")
        add_message("room_1", "info", "second")

        output = read_output("room_1", since=1)

        self.assertEqual([{"level": "info", "text": "second"}], output["messages"])
        self.assertEqual(2, output["next"])

    def test_the_same_message_twice_in_a_row_is_kept_once(self):
        start_command("room_1")
        add_message("room_1", "error", "Unknown experiment")
        add_message("room_1", "error", "Unknown experiment")

        self.assertEqual(1, len(read_output("room_1", since=0)["messages"]))

    def test_a_command_without_errors_is_completed(self):
        start_command("room_1")
        add_message("room_1", "warning", "No files found.")

        finish_command("room_1")

        output = read_output("room_1", since=0)
        self.assertEqual("completed", output["status"])
        self.assertEqual(
            {"level": "info", "text": "Command completed."}, output["messages"][-1]
        )

    def test_a_command_with_an_error_failed(self):
        start_command("room_1")
        add_message("room_1", "error", "Unknown experiment")

        finish_command("room_1")

        output = read_output("room_1", since=0)
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            {"level": "error", "text": "Command failed."}, output["messages"][-1]
        )

    def test_a_new_command_starts_without_the_old_messages(self):
        start_command("room_1")
        add_message("room_1", "error", "old")
        finish_command("room_1")

        start_command("room_1")

        self.assertEqual(
            {"messages": [], "next": 0, "status": "running"},
            read_output("room_1", since=0),
        )

    def test_before_the_command_starts_there_is_no_status(self):
        self.assertEqual(
            {"messages": [], "next": 0, "status": None},
            read_output("room_1", since=0),
        )

    def test_without_room_name_nothing_is_stored(self):
        start_command(None)
        add_message(None, "error", "Only in the log")
        finish_command(None)

        self.assertEqual(
            {"messages": [], "next": 0, "status": None},
            read_output("None", since=0),
        )

    def test_a_command_with_very_many_messages_stops_collecting_them(self):
        start_command("room_1")
        for number in range(MAX_MESSAGES + 10):
            add_message("room_1", "info", f"line {number}")

        finish_command("room_1")

        messages = read_output("room_1", since=0)["messages"]
        self.assertEqual(MAX_MESSAGES + 2, len(messages))
        self.assertEqual(
            {"level": "warning", "text": TOO_MANY_MESSAGES}, messages[MAX_MESSAGES]
        )
        self.assertEqual({"level": "info", "text": "Command completed."}, messages[-1])


class InterruptedCommandTest(SimpleTestCase):
    def setUp(self):
        cache.clear()

    def test_a_command_that_was_running_when_the_worker_restarted_failed(self):
        start_command("room_1")
        register_running_command("room_1", "celery@worker_1")
        add_message("room_1", "info", "Processing file a.csv...")

        fail_interrupted_commands("celery@worker_1")

        output = read_output("room_1", since=0)
        self.assertEqual("failed", output["status"])
        self.assertEqual(
            [
                {"level": "info", "text": "Processing file a.csv..."},
                {"level": "error", "text": INTERRUPTED_MESSAGE},
                {"level": "error", "text": "Command failed."},
            ],
            output["messages"],
        )

    def test_a_finished_command_is_not_touched(self):
        start_command("room_1")
        register_running_command("room_1", "celery@worker_1")
        finish_command("room_1")

        fail_interrupted_commands("celery@worker_1")

        output = read_output("room_1", since=0)
        self.assertEqual("completed", output["status"])
        self.assertEqual(1, len(output["messages"]))

    def test_a_command_waiting_in_the_queue_is_not_touched(self):
        # The view has set "running", but no worker has started the command yet
        start_command("room_1")

        fail_interrupted_commands("celery@worker_1")

        self.assertEqual(
            {"messages": [], "next": 0, "status": "running"},
            read_output("room_1", since=0),
        )

    def test_a_command_of_another_worker_is_not_touched(self):
        start_command("room_1")
        register_running_command("room_1", "celery@worker_2")

        fail_interrupted_commands("celery@worker_1")

        self.assertEqual("running", read_output("room_1", since=0)["status"])


class ErrorTextTest(SimpleTestCase):
    def test_a_command_error_is_shown_as_it_is(self):
        self.assertEqual(
            "Project P1 does not exist.",
            error_text(CommandError("Project P1 does not exist.")),
        )

    def test_any_other_error_shows_its_type(self):
        self.assertEqual("KeyError: 'NAME'", error_text(KeyError("NAME")))
        # The text alone would be empty
        self.assertEqual("AssertionError", error_text(AssertionError()))
