"""
PDF reports, report files and measurement CSV downloads.
"""

import os
import subprocess
from django.conf import settings
from django.http import HttpResponse
from rest_framework import status
from django.http import JsonResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from helpers.paths import path_inside
from helpers.logger import logger
from ldm.ldm import get_experiment_measurements


def notebook_path_inside_reports(path: str) -> str:
    """A notebook or report path, if it is inside settings.NOTEBOOKS_ROOT."""
    return path_inside(path, settings.NOTEBOOKS_ROOT)


# These views read and write report files, so they are only for logged in users
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def generate_pdf_report(request):
    """Runs the notebook that writes the PDF report of an experiment."""
    notebook_path = notebook_path_inside_reports(
        request.data.get("notebook_path")
        or f"{settings.NOTEBOOKS_ROOT}/input/general.ipynb"
    )
    command = [
        "python",
        "/root/.ipython/profile_default/report_generator.py",
        "--notebook_path",
        notebook_path,
        "--experiment",
        request.data.get("experiment"),
        "--label",
        request.data.get("label"),
        "--pos",
        request.data.get("selected_pos"),
        "--neg",
        request.data.get("selected_neg"),
    ]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as error:
        logger.exception("The report notebook failed")
        return JsonResponse({"error": str(error)}, status=500)
    return JsonResponse({"status": "Report generated successfully"})


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def list_files(request):
    """The report files of an experiment, e.g. its PDF reports."""
    experiment = request.data.get("experiment")
    notebooks_dir = notebook_path_inside_reports(
        request.data.get("notebooks_dir")
        or f"{settings.NOTEBOOKS_ROOT}/output/{experiment}"
    )
    file_format = request.data.get("file_format") or ".pdf"

    notebooks = []
    if os.path.exists(notebooks_dir):
        notebooks = [
            os.path.join(notebooks_dir, file)
            for file in os.listdir(notebooks_dir)
            if file.endswith(file_format)
        ]
    return Response({"notebooks": notebooks}, status=status.HTTP_200_OK)


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def download_pdf_report(request):
    """One report file, for the download in the browser."""
    path = notebook_path_inside_reports(request.data.get("path"))
    if not os.path.exists(path):
        return JsonResponse({"error": "File not found"}, status=404)

    response = HttpResponse(open(path, "rb"), content_type="application/pdf")
    response["Content-Disposition"] = f'attachment; filename="{os.path.basename(path)}"'
    return response


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def download_csv_data(request):
    """The measurements of an experiment as a csv file."""
    label = request.data.get("label")
    experiment = request.data.get("experiment")
    measurement_type = request.data.get("type")

    data_frame = get_experiment_measurements(
        experiment, label, measurement_type, csv=True
    )
    response = HttpResponse(content_type="text/csv")
    suffix = "_meas" if measurement_type == "main" else "_comp"
    response["Content-Disposition"] = f'attachment; filename="{label}_{suffix}.csv"'
    data_frame.to_csv(response, index=False)
    return response
