import os.path
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json
from io import StringIO
import sys
from django.core import management
from django.core.cache import cache

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
