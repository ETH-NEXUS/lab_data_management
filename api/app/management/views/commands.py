"""
Running a command from the management page and reading its output.

The views of the management page are for logged in users only. As DRF views
they also check the CSRF token of POST requests, which the UI sends. DRF reads
the request body for that check, so the views use request.data (or request.POST
for an upload), never request.body.
"""

import re

from django.http import JsonResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.exceptions import ValidationError
from rest_framework.permissions import IsAuthenticated

from importer.command_output import read_output, start_command
from importer.tasks import run_management_command
from management.paths import check_command_paths

# A room name is made by the page, e.g. "12_1726563600000"
ROOM_NAME = re.compile(r"\w{1,64}")


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def run_command(request):
    form_data = request.data.get("form_data")
    if not isinstance(form_data, dict):
        raise ValidationError("The request has no form_data.")
    check_command_paths(form_data)

    room_name = form_data.get("room_name")
    if not isinstance(room_name, str) or not ROOM_NAME.fullmatch(room_name):
        raise ValidationError(f"This is not a room name: {room_name}")
    start_command(room_name)
    # The command runs in the celery container; the page reads its output
    # through long_polling while it runs
    run_management_command.delay(form_data)

    return JsonResponse({"status": "ok"})


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def long_polling(request, room_name):
    """
    The new messages of a running command, from position `since` on.
    Example: GET /api/long_polling/12_1726/?since=3
    -> {"messages": [{"level": "error", "text": "..."}], "next": 4, "status": "failed"}
    """
    since = request.GET.get("since", "0")
    since = int(since) if since.isdigit() else 0
    return JsonResponse(read_output(room_name, since))
