from django.http import JsonResponse, HttpResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.exceptions import ValidationError
from rest_framework.permissions import IsAuthenticated
import os
from django.http import Http404
from django.conf import settings

from importer.command_output import read_output, start_command
from management.paths import check_command_paths, data_path
from importer.tasks import run_management_command
from chardet.universaldetector import UniversalDetector
from contextlib import redirect_stderr


def list_files(start_path):
    def walk(path):
        data = {
            "type": "directory",
            "name": os.path.basename(path),
            "children": [],
            "path": path,
        }

        if not ".snapshots" in path:
            for entry in os.scandir(path):
                if entry.is_file():
                    data["children"].append(
                        {"type": "file", "name": entry.name, "path": entry.path}
                    )
                elif entry.is_dir():
                    data["children"].append(walk(entry.path))
        return data

    return walk(start_path)


# The views of the management page are for logged in users only. As DRF views
# they also check the CSRF token of POST requests, which the UI sends. DRF reads
# the request body for that check, so the views use request.data, not request.body.
@api_view(["GET"])
@permission_classes([IsAuthenticated])
def directory_content(request):
    content = list_files(settings.MANAGEMENT_DATA_ROOT)
    return JsonResponse({"directory_content": content})


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def run_command(request):
    form_data = request.data.get("form_data")
    if not isinstance(form_data, dict):
        raise ValidationError("The request has no form_data.")
    check_command_paths(form_data)

    room_name = form_data.get("room_name")
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


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def delete_file(request):
    path = data_path(request.data.get("path"))
    if os.path.exists(path):
        os.remove(path)
    return JsonResponse({"status": "ok"})


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def download_file(request):
    file_path = data_path(request.data.get("file_path"))
    if not os.path.exists(file_path):
        raise Http404("File not found")

    response = HttpResponse(
        open(file_path, "rb"), content_type="application/octet-stream"
    )
    response["Content-Disposition"] = (
        f'attachment; filename="{os.path.basename(file_path)}"'
    )
    return response


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def upload_file(request):
    directory_path = request.POST.get("directory_path")
    uploaded_file = request.FILES.get("file")
    if not (directory_path and uploaded_file):
        raise ValidationError("The request has no directory path or no file.")

    directory_path = data_path(directory_path)
    os.makedirs(directory_path, exist_ok=True)
    # Only the name of the uploaded file, never a path it may carry
    file_path = os.path.join(directory_path, os.path.basename(uploaded_file.name))
    with open(file_path, "wb+") as destination:
        for chunk in uploaded_file.chunks():
            destination.write(chunk)

    return JsonResponse(
        {"message": "File uploaded successfully", "file_path": file_path}
    )


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def get_file_content(request):
    file_path = data_path(request.data.get("file_path"))
    if not os.path.exists(file_path):
        raise Http404("File not found")

    with redirect_stderr(None):
        detector = UniversalDetector()
        with open(file_path, "rb") as file:
            for line in file:
                detector.feed(line)
                if detector.done:
                    break
            detector.close()
        encoding = detector.result.get("encoding")

    with open(file_path, "r", encoding=encoding) as file:
        return JsonResponse({"content": file.read()})
