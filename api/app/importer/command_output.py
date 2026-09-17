"""
The output of a command that runs from the management page.

While the command runs, its messages are kept in the cache (Redis), together with
their level. The page asks `long_polling` for the new messages and the status.

Stored messages example:
[{"level": "info", "text": "Processing file /data/a.csv..."},
 {"level": "error", "text": "Cannot read the measurement date of /data/b.txt ..."}]
"""

from django.core.cache import cache

# The output of a command is kept this long; Redis removes it afterwards.
OUTPUT_TIMEOUT_SECONDS = 24 * 60 * 60

RUNNING = "running"
COMPLETED = "completed"
FAILED = "failed"


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
