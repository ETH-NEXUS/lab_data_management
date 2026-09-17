import os.path
from django.http import JsonResponse, HttpResponse
from django.views.decorators.csrf import csrf_exempt
import json
import os
from django.http import Http404
from importer.command_output import read_output, start_command
from importer.tasks import run_management_command
from chardet.universaldetector import UniversalDetector
from contextlib import redirect_stderr


def list_files(start_path):
    def walk(path, children=None):
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
                    data["children"].append(walk(entry.path, children))
        return data

    return walk(start_path)


def directory_content(request, start_path="/data"):
    content = list_files(start_path)
    return JsonResponse({"directory_content": content})


@csrf_exempt
def run_command(request):
    if request.method == "POST":
        body_unicode = request.body.decode("utf-8")
        body_data = json.loads(body_unicode)
        form_data = body_data.get("form_data")

        room_name = form_data.get("room_name")
        start_command(room_name)
        # The command runs in the celery container; the page reads its output
        # through long_polling while it runs
        run_management_command.delay(form_data)

    return JsonResponse({"status": "ok"})


def long_polling(request, room_name):
    """
    The new messages of a running command, from position `since` on.
    Example: GET /api/long_polling/12_1726/?since=3
    -> {"messages": [{"level": "error", "text": "..."}], "next": 4, "status": "failed"}
    """
    since = request.GET.get("since", "0")
    since = int(since) if since.isdigit() else 0
    return JsonResponse(read_output(room_name, since))


@csrf_exempt
def delete_file(request):
    if request.method == "POST":
        body_unicode = request.body.decode("utf-8")
        body_data = json.loads(body_unicode)
        path = body_data.get("path")
        if os.path.exists(path):
            os.remove(path)
        return JsonResponse({"status": "ok"})

    return JsonResponse({"status": "error"})


@csrf_exempt
def download_file(request):
    if request.method == "POST":
        body_unicode = request.body.decode("utf-8")
        body_data = json.loads(body_unicode)
        file_path = body_data.get("file_path")
        print("file_path", file_path)

        if not file_path:
            raise Http404("File path not provided")

        if os.path.exists(file_path):
            try:
                file = open(file_path, "rb")
            except IOError:
                raise Http404("File not found")

            response = HttpResponse(file, content_type="application/octet-stream")
            response[
                "Content-Disposition"
            ] = f'attachment; filename="{os.path.basename(file_path)}"'
            return response
        else:
            raise Http404("File not found")


@csrf_exempt
def upload_file(request):
    if request.method == "POST":
        directory_path = request.POST.get("directory_path")
        uploaded_file = request.FILES.get("file")

        if not (directory_path and uploaded_file):
            return JsonResponse(
                {"status": "error", "error": "Directory path or file not provided"}
            )
        os.makedirs(directory_path, exist_ok=True)
        file_path = os.path.join(directory_path, uploaded_file.name)

        with open(file_path, "wb+") as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        return JsonResponse(
            {"message": "File uploaded successfully", "file_path": file_path}
        )


@csrf_exempt
def get_file_content(request):
    if request.method == "POST":
        body_unicode = request.body.decode("utf-8")
        body_data = json.loads(body_unicode)
        file_path = body_data.get("file_path")

        if os.path.exists(file_path):
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
                content = file.read()
            return JsonResponse({"content": content})
        else:
            raise Http404("File not found")
    else:
        return JsonResponse({"error": "Invalid request method"})
