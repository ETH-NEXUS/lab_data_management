"""
The output of a command that runs from the management page.

While the command runs, its messages are kept in the cache (Redis), together with
their level. The page asks `long_polling` for the new messages and the status.

Stored messages example:
[{"level": "info", "text": "Processing file /data/a.csv..."},
 {"level": "error", "text": "/data/b.txt was not mapped, nothing of it was stored: ..."}]
"""

from django.core.cache import cache
from django.core.management.base import CommandError

# The output of a command is kept this long; Redis removes it afterwards.
OUTPUT_TIMEOUT_SECONDS = 24 * 60 * 60

RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"

# The rooms of the commands that a worker is running right now, e.g. ["12_1726563600000"]
RUNNING_COMMANDS_KEY = "running_commands"
INTERRUPTED_MESSAGE = (
    "The command was interrupted, because the command worker was restarted. "
    "Check what was stored before you run it again."
)


def messages_key(room_name: str) -> str:
    return f"command_messages_{room_name}"


def status_key(room_name: str) -> str:
    return f"command_status_{room_name}"


def start_command(room_name: str | None) -> None:
    """A new command starts: no messages yet, status "running"."""
    if not room_name:
        return
    cache.set(messages_key(room_name), [], OUTPUT_TIMEOUT_SECONDS)
    cache.set(status_key(room_name), RUNNING, OUTPUT_TIMEOUT_SECONDS)


def register_running_command(room_name: str | None) -> None:
    """The worker starts the command. It is listed until finish_command."""
    if not room_name:
        return
    running = cache.get(RUNNING_COMMANDS_KEY, [])
    if room_name not in running:
        running.append(room_name)
    cache.set(RUNNING_COMMANDS_KEY, running, OUTPUT_TIMEOUT_SECONDS)


def fail_interrupted_commands() -> None:
    """
    Called when the worker starts: a command that is still listed as running
    was stopped by the restart, so it gets an error and the status "failed".
    """
    for room_name in cache.get(RUNNING_COMMANDS_KEY, []):
        add_message(room_name, "error", INTERRUPTED_MESSAGE)
        finish_command(room_name)


def add_message(room_name: str | None, level: str, text: str) -> None:
    """Adds one message, e.g. add_message("12_1726", "error", "File not found")."""
    if not room_name:
        return
    messages = cache.get(messages_key(room_name), [])
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
    has_errors = any(message["level"] == "error" for message in messages)
    if has_errors:
        add_message(room_name, "error", "Command failed.")
        status = FAILED
    else:
        add_message(room_name, "info", "Command completed.")
        status = COMPLETED
    cache.set(status_key(room_name), status, OUTPUT_TIMEOUT_SECONDS)

    running = cache.get(RUNNING_COMMANDS_KEY, [])
    if room_name in running:
        running.remove(room_name)
        cache.set(RUNNING_COMMANDS_KEY, running, OUTPUT_TIMEOUT_SECONDS)


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
