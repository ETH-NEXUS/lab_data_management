"""
Starting the statistical analysis of an experiment, and downloading its results.
The output of a running analysis is read through long_polling, like a command
of the management page.
"""

import os
import re
import time

from django.http import FileResponse, JsonResponse
from rest_framework.decorators import api_view, permission_classes
from rest_framework.exceptions import NotFound, ValidationError
from rest_framework.permissions import IsAuthenticated

from analysis.tasks import ANALYSIS_FOLDER, run_analysis
from core.models import PlateInfo
from importer.command_output import read_output, start_command
from management.views.commands import check_room_name

# A result is the zip of one run, e.g. "20260929-125825_single_Lum_CTG.zip"
RESULT_NAME = re.compile(r"[\w-]+\.zip")


@api_view(["POST"])
@permission_classes([IsAuthenticated])
def start_analysis(request):
    """
    Starts the analysis in the celery container.

    Accepted data example:
    {"experiment_id": 105, "label": "Lum_CTG", "analysis_type": "single",
     "settings": {"condi_yes": "irradiated"}, "room_name": "3_1727600000000"}
    """
    form_data = dict(request.data)
    room_name = form_data.get("room_name")
    check_room_name(room_name)
    # A second start with the same room would clear the output of the running one
    if read_output(room_name, 0)["status"] is not None:
        raise ValidationError(f"The analysis {room_name} was started already.")
    # The worker does not start an analysis that waited too long for it
    form_data["queued_at"] = time.time()
    start_command(room_name)
    run_analysis.delay(form_data)
    return JsonResponse({"status": "ok"})


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def list_conditions(request):
    """
    The conditions in the saved plate information of one measurement, for the
    selectivity analysis.
    Example: GET /api/analysis/conditions/?experiment_id=86&label=Lum
    -> {"conditions": ["irradiated", "not irradiated"]}
    """
    conditions = (
        PlateInfo.objects.filter(
            experiment_id=experiment_id(request), label=request.GET.get("label", "")
        )
        .exclude(condition="")
        .values_list("condition", flat=True)
        .distinct()
    )
    return JsonResponse({"conditions": sorted(conditions)})


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def list_analysis_results(request):
    """
    The result zips of an experiment, newest first.
    Example: GET /api/analysis/results/?experiment_id=105
    -> {"results": ["20260929-125825_single_Lum_CTG.zip"]}
    """
    experiment_folder = os.path.join(ANALYSIS_FOLDER, experiment_id(request))
    results = []
    if os.path.isdir(experiment_folder):
        # The names start with the date and time of the run, so this sorts by time
        for file_name in sorted(os.listdir(experiment_folder), reverse=True):
            if RESULT_NAME.fullmatch(file_name):
                results.append(file_name)
    return JsonResponse({"results": results})


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def download_analysis_result(request):
    """
    One result zip.
    Example: GET /api/analysis/download/?experiment_id=105&name=20260929-125825_single_Lum_CTG.zip
    """
    name = request.GET.get("name", "")
    if not RESULT_NAME.fullmatch(name):
        raise ValidationError(f"This is not the name of a result: {name}")
    zip_path = os.path.join(ANALYSIS_FOLDER, experiment_id(request), name)
    if not os.path.exists(zip_path):
        raise NotFound(f"The result {name} does not exist.")
    return FileResponse(open(zip_path, "rb"), as_attachment=True, filename=name)


def experiment_id(request) -> str:
    """The experiment id of the request, e.g. "105"."""
    value = request.GET.get("experiment_id", "")
    if not value.isdigit():
        raise ValidationError(f"This is not an experiment id: {value}")
    return value
