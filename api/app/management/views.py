import os.path
from django.http import JsonResponse, HttpResponse
from django.views.decorators.csrf import csrf_exempt
import json
import os
from django.conf import settings
from django.http import FileResponse, Http404
from django.core import management
from django.core.cache import cache
from django.core.files.storage import default_storage

from importer.helper import message


def list_files(start_path):
    def walk(path, children=None):
        data = {
            "type": "directory",
            "name": os.path.basename(path),
            "children": [],
            "path": path,
        }
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
    # output = StringIO()
    # sys.stdout = output

    if request.method == "POST":
        body_unicode = request.body.decode("utf-8")
        body_data = json.loads(body_unicode)
        form_data = body_data.get("form_data")

        if form_data.get("command") == "map":
            machine = form_data.get("machine")
            if machine in ["echo", "m1000"]:
                kwargs = {
                    "path": form_data.get("path"),
                    "mapping_file": form_data.get("mapping_file"),
                    "debug": False,
                    "create_missing_plates": True,
                    "experiment_name": form_data.get("experiment_name"),
                    "room_name": form_data.get("room_name"),
                }
                management.call_command("map", machine, **kwargs)
        elif form_data.get("command") == "import":
            what = form_data.get("what")
            kwargs = {
                "mapping_file": form_data.get("mapping_file"),
                "input_file": form_data.get("input_file"),
                "debug": False,
                "library_name": form_data.get("library_name")
                if form_data.get("library_name")
                else None,
                "template_name": form_data.get("template_name")
                if form_data.get("template_name")
                else None,
                "room_name": form_data.get("room_name"),
            }
            management.call_command("import", what, **kwargs)
        message("Command completed.", "info", form_data.get("room_name"))
        cache.set(f"command_status_{form_data.get('room_name')}", "completed")

    # sys.stdout = sys.__stdout__
    # output_str = output.getvalue()
    # output.close()

    return JsonResponse({"status": "ok"})


def long_polling(request, room_name):
    output_key = f"command_output_{room_name}"
    status_key = f"command_status_{room_name}"
    output = cache.get(output_key)
    status = cache.get(status_key)

    if output:
        cache.delete(output_key)
        return JsonResponse({"message": output, "status": status})
    else:
        return JsonResponse({"message": None, "status": status})


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
