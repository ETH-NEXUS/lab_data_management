import os.path
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json
from io import StringIO
import sys
from django.core import management


def list_files(start_path):
    def walk(path, children=None):
        data = {"type": "directory", "name": os.path.basename(path), "children": []}
        for entry in os.scandir(path):
            if entry.is_file():
                data["children"].append({"type": "file", "name": entry.name})
            elif entry.is_dir():
                data["children"].append(walk(entry.path, children))
        return data

    return walk(start_path)

def directory_content(request, start_path="/data"):
    content = list_files(start_path)
    return JsonResponse({"directory_content": content})


@csrf_exempt
def run_command(request):
    output = StringIO()
    sys.stdout = output

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
                    "debug": form_data.get("debug"),
                    "create_missing_plates": True,
                    "experiment_name": form_data.get("experiment_name"),
                    "stdout": output,
                }
                management.call_command("map", machine, **kwargs)

    sys.stdout = sys.__stdout__
    output_str = output.getvalue()
    output.close()

    return JsonResponse({"command_output": output_str})

