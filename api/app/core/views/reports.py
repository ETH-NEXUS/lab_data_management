"""
PDF reports, report files and measurement CSV downloads.
"""

import json
import os
import subprocess
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
from rest_framework import status
from django.http import JsonResponse
from rest_framework.decorators import api_view
from rest_framework.response import Response
from ldm.ldm import get_experiment_measurements


@csrf_exempt
def generate_pdf_report(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            notebook_path = data.get("notebook_path")
            experiment = data.get("experiment")
            label = data.get("label")
            selected_pos = data.get("selected_pos")
            selected_neg = data.get("selected_neg")
            if not notebook_path:

                notebook_path = "/notebooks/input/general.ipynb"

            cmd = [
                "python",
                "/root/.ipython/profile_default/report_generator.py",
                "--notebook_path",
                notebook_path,
                "--experiment",
                experiment,
                "--label",
                label,
                "--pos",
                selected_pos,
                "--neg",
                selected_neg,
            ]
            subprocess.run(cmd, check=True)
            return JsonResponse({"status": "Report generated successfully"})
        else:
            return JsonResponse({"error": "Invalid request method"}, status=400)
    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
@api_view(["POST"])
def list_files(request):
    try:
        experiment = request.data.get("experiment")
        notebooks_dir = request.data.get("notebooks_dir")
        file_format = request.data.get("file_format")
        print(experiment, notebooks_dir, file_format)
        if not notebooks_dir:
            notebooks_dir = f"/notebooks/output/{experiment}"
        if not file_format:
            file_format = ".pdf"
        notebooks = []
        if os.path.exists(notebooks_dir):
            for file in os.listdir(notebooks_dir):
                if file.endswith(file_format):
                    full_path = os.path.join(notebooks_dir, file)
                    notebooks.append(full_path)
        return Response({"notebooks": notebooks}, status=status.HTTP_200_OK)
    except Exception as e:
        print(e)
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


@csrf_exempt
def download_pdf_report(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            path = data.get("path")

            if not path:
                return JsonResponse({"error": "Path not provided"}, status=400)
            if not os.path.exists(path):
                return JsonResponse({"error": "File not found"}, status=404)
            with open(path, "rb") as f:
                response = HttpResponse(f, content_type="application/pdf")
                response[
                    "Content-Disposition"
                ] = f'attachment; filename="{os.path.basename(path)}"'
                return response
        else:
            return JsonResponse({"error": "Invalid request method"}, status=400)
    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)


@csrf_exempt
def download_csv_data(request):
    try:
        if request.method == "POST":
            data = json.loads(request.body.decode("utf-8"))
            label = data.get("label")
            experiment = data.get("experiment")
            type = data.get("type")
            df = get_experiment_measurements(experiment, label, type, csv=True)
            response = HttpResponse(content_type="text/csv")
            suffix = "_meas" if type == "main" else "_comp"
            response[
                "Content-Disposition"
            ] = f'attachment; filename="{label}_{suffix}.csv"'
            df.to_csv(response, index=False)
            return response
        return None

    except Exception as e:
        print("EXCEPTION")
        print(e)
        return JsonResponse({"error": str(e)}, status=500)
