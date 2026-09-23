"""
The files of the data folder on the management page: the folder tree, reading,
downloading, uploading and deleting a file.

The views of the management page are for logged in users only. As DRF views
they also check the CSRF token of POST requests, which the UI sends. DRF reads
the request body for that check, so the views use request.data (or request.POST
for an upload), never request.body.
"""

import os

from django.conf import settings
from django.http import FileResponse, Http404, JsonResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.exceptions import ValidationError
from rest_framework.permissions import IsAuthenticated

from importer.mappers.base import detect_encoding
from management.paths import data_path


def folder_tree(path: str) -> dict:
    """
    A folder with everything in it, e.g.
    {"type": "directory", "name": "echo", "path": "/data/echo", "children": [
        {"type": "file", "name": "run.xml", "path": "/data/echo/run.xml"}]}
    """
    tree = {
        "type": "directory",
        "name": os.path.basename(path),
        "children": [],
        "path": path,
    }
    # The snapshots of the lab shares hold old copies of every file; they are
    # listed as a folder, but not opened
    if ".snapshots" in path:
        return tree

    for entry in os.scandir(path):
        if entry.is_file():
            tree["children"].append(
                {"type": "file", "name": entry.name, "path": entry.path}
            )
        # A link to another folder is not followed: the lab shares
        # contain links that would send this into a circle
        elif entry.is_dir(follow_symlinks=False):
            tree["children"].append(folder_tree(entry.path))
    return tree


def existing_file(path: str) -> str:
    """The path of a file inside the data folder that can be read."""
    path = data_path(path)
    if os.path.isdir(path):
        raise ValidationError(f"{path} is a folder, not a file.")
    if not os.path.exists(path):
        raise Http404("File not found")
    return path


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def directory_content(request):
    content = folder_tree(settings.MANAGEMENT_DATA_ROOT)
    return JsonResponse({"directory_content": content})


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def delete_file(request):
    path = data_path(request.data.get("path"))
    if os.path.isdir(path):
        raise ValidationError(f"{path} is a folder, only files can be deleted.")
    if os.path.exists(path):
        os.remove(path)
    return JsonResponse({"status": "ok"})


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def download_file(request):
    file_path = existing_file(request.data.get("file_path"))
    # FileResponse sends the file in parts and closes it afterwards
    return FileResponse(
        open(file_path, "rb"),
        as_attachment=True,
        filename=os.path.basename(file_path),
        content_type="application/octet-stream",
    )


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
    file_name = os.path.basename(uploaded_file.name)
    file_path = os.path.join(directory_path, file_name)
    # "x" creates the file only if there is none yet, so an upload never
    # replaces a file that is already there (or that is being written right now)
    try:
        destination = open(file_path, "xb")
    except FileExistsError:
        raise ValidationError(
            f"A file named {file_name} already exists in {directory_path}. "
            "Delete it first or rename the new file."
        )
    with destination:
        for chunk in uploaded_file.chunks():
            destination.write(chunk)

    return JsonResponse(
        {"message": "File uploaded successfully", "file_path": file_path}
    )


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def get_file_content(request):
    file_path = existing_file(request.data.get("file_path"))
    encoding = detect_encoding(file_path)
    with open(file_path, "r", encoding=encoding) as file:
        return JsonResponse({"content": file.read()})
