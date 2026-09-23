"""
The output of a command that runs from the management page.

While the command runs, its messages are kept in the cache (Redis), together with
their level. The page asks `long_polling` for the new messages and the status.

Stored messages example:
[{"level": "info", "text": "Processing file /data/a.csv..."},
 {"level": "error", "text": "/data/b.txt was not mapped, nothing of it was stored: ..."}]
"""

import time
from contextlib import contextmanager

from django.core.cache import cache
from django.core.management.base import CommandError

# The output of a command is kept this long; Redis removes it afterwards.
OUTPUT_TIMEOUT_SECONDS = 24 * 60 * 60

RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"

# A command does not show more messages than this; the rest is only in the log
MAX_MESSAGES = 1000
TOO_MANY_MESSAGES = (
    f"More than {MAX_MESSAGES} messages: the next ones are only in the server log."
)

# The worker of every command that is running right now, by room:
# {"12_1726563600000": "celery@celery"}
RUNNING_COMMANDS_KEY = "running_commands"
RUNNING_COMMANDS_LOCK_KEY = "running_commands_lock"
LOCK_TIMEOUT_SECONDS = 10
INTERRUPTED_MESSAGE = (
    "The command was interrupted, because the command worker was restarted. "
    "Check what was stored before you run it again."
)
LOST_PROCESS_MESSAGE = (
    "The command was stopped, because the process that ran it was killed "
    "(for example, it used too much memory). The file it was reading was not "
    "stored; check what was stored before you run it again."
)


def messages_key(room_name: str) -> str:
    return f"command_messages_{room_name}"


def error_flag_key(room_name: str) -> str:
    return f"command_had_error_{room_name}"


def status_key(room_name: str) -> str:
    return f"command_status_{room_name}"


def start_command(room_name: str | None) -> None:
    """A new command starts: no messages yet, status "running"."""
    if not room_name:
        return
    cache.set(messages_key(room_name), [], OUTPUT_TIMEOUT_SECONDS)
    cache.set(error_flag_key(room_name), False, OUTPUT_TIMEOUT_SECONDS)
    cache.set(status_key(room_name), RUNNING, OUTPUT_TIMEOUT_SECONDS)


@contextmanager
def running_commands_lock():
    """
    Holds the list of running commands while it is read and written again, so
    that two commands starting at the same moment do not overwrite each other.
    """
    while not cache.add(RUNNING_COMMANDS_LOCK_KEY, "1", LOCK_TIMEOUT_SECONDS):
        time.sleep(0.05)
    try:
        yield
    finally:
        cache.delete(RUNNING_COMMANDS_LOCK_KEY)


def register_running_command(room_name: str | None, worker_name: str) -> None:
    """The worker starts the command. It is listed until finish_command."""
    if not room_name:
        return
    with running_commands_lock():
        running = cache.get(RUNNING_COMMANDS_KEY, {})
        running[room_name] = worker_name
        cache.set(RUNNING_COMMANDS_KEY, running, OUTPUT_TIMEOUT_SECONDS)


def forget_running_command(room_name: str) -> None:
    """The command has ended, so its worker is not running it anymore."""
    with running_commands_lock():
        running = cache.get(RUNNING_COMMANDS_KEY, {})
        if running.pop(room_name, None) is not None:
            cache.set(RUNNING_COMMANDS_KEY, running, OUTPUT_TIMEOUT_SECONDS)


def fail_interrupted_commands(worker_name: str) -> None:
    """
    Called when a worker starts: a command that this worker was running was
    stopped by the restart, so it gets an error and the status "failed".
    Commands of another worker are left alone.
    """
    running = cache.get(RUNNING_COMMANDS_KEY, {})
    for room_name in [
        room for room, worker in running.items() if worker == worker_name
    ]:
        add_message(room_name, "error", INTERRUPTED_MESSAGE)
        finish_command(room_name)


def fail_lost_command(room_name: str | None) -> None:
    """
    The process that ran the command was killed in the middle of it, so the
    command could not end itself: it gets an error and the status "failed".
    """
    add_message(room_name, "error", LOST_PROCESS_MESSAGE)
    finish_command(room_name)


def add_message(room_name: str | None, level: str, text: str) -> None:
    """Adds one message, e.g. add_message("12_1726", "error", "File not found")."""
    if not room_name:
        return
    if level == "error":
        # Remembered separately: after the limit below, the level is not kept
        cache.set(error_flag_key(room_name), True, OUTPUT_TIMEOUT_SECONDS)

    messages = cache.get(messages_key(room_name), [])
    if len(messages) >= MAX_MESSAGES:
        # The same message twice in a row is dropped below, so this line is added once
        level, text = "warning", TOO_MANY_MESSAGES
    new_message = {"level": level, "text": text}
    # A command that reports an error and then raises it would show it twice
    if messages and messages[-1] == new_message:
        return
    messages.append(new_message)
    cache.set(messages_key(room_name), messages, OUTPUT_TIMEOUT_SECONDS)


def finish_command(room_name: str | None) -> None:
    """
    The command has ended. It "failed" if at least one message is an error,
    otherwise it is "completed". A last message says which one.
    """
    if not room_name:
        return
    messages = cache.get(messages_key(room_name), [])
    has_errors = cache.get(error_flag_key(room_name)) or any(
        message["level"] == "error" for message in messages
    )
    if has_errors:
        messages.append({"level": "error", "text": "Command failed."})
        status = FAILED
    else:
        messages.append({"level": "info", "text": "Command completed."})
        status = COMPLETED
    cache.set(messages_key(room_name), messages, OUTPUT_TIMEOUT_SECONDS)
    cache.set(status_key(room_name), status, OUTPUT_TIMEOUT_SECONDS)
    forget_running_command(room_name)


def read_output(room_name: str, since: int) -> dict:
    """
    The messages from position `since` on, the position to ask from next time,
    and the status (None while the command has not started yet).

    Returned data example:
    {"messages": [{"level": "info", "text": "Command completed."}],
     "next": 5, "status": "completed"}
    """
    messages = cache.get(messages_key(room_name), [])
    return {
        "messages": messages[since:],
        "next": len(messages),
        "status": cache.get(status_key(room_name)),
    }


def error_text(error: Exception) -> str:
    """
    The text of an error for the management page. A CommandError is written
    for people, so only its text is shown. Any other error is unexpected: its
    type is shown too, because its text alone can be empty or unclear.

    CommandError("Project P1 does not exist.") -> "Project P1 does not exist."
    KeyError("NAME") -> "KeyError: 'NAME'"
    AssertionError() -> "AssertionError"
    """
    if isinstance(error, CommandError):
        return str(error)
    error_type = type(error).__name__
    if not str(error):
        return error_type
    return f"{error_type}: {error}"
