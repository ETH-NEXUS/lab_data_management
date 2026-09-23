"""
Projects of this app can be linked to a project in Harvest, the time tracking
of the lab. They then take over its name and notes.
"""

import logging

import requests
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import get_object_or_404
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated

from core.models import Project

from .harvest_client import HarvestClient

logger = logging.getLogger(__name__)

client = HarvestClient(settings.HARVEST_ACCESS_TOKEN, settings.HARVEST_ACCOUNT_ID)

HARVEST_UNAVAILABLE = "Harvest could not be asked, please try again later."
HARVEST_NOT_SET_UP = "Harvest is not set up on this server."


def harvest_unavailable(error: requests.RequestException) -> JsonResponse:
    """The answer when Harvest cannot be reached or answers with an error."""
    logger.warning("Harvest request failed: %s", error)
    return JsonResponse({"error": HARVEST_UNAVAILABLE}, status=502)


# Both views are for logged in users only. The update stays a GET request,
# because the UI sends it as one.
@api_view(["GET"])
@permission_classes([IsAuthenticated])
def harvest_projects(request):
    """
    The Harvest projects the user can choose from, e.g.
    {"projects": [{"id": 7, "name": "Screening 2026", "notes": "...", ...}]}
    """
    if settings.HARVEST_ACCESS_TOKEN is None:
        return JsonResponse({"projects": []})

    try:
        projects = client.get("projects")["projects"]
    except requests.RequestException as error:
        return harvest_unavailable(error)

    name_filter = settings.HARVEST_PROJECT_FILTER
    if name_filter:
        projects = [project for project in projects if name_filter in project["name"]]
    return JsonResponse({"projects": projects})


@api_view(["GET"])
@permission_classes([IsAuthenticated])
def update_harvest_info(request, project_id):
    """Takes over the name and the notes of the linked Harvest project."""
    project = get_object_or_404(Project, id=project_id)
    if not project.harvest_id:
        return JsonResponse({"success": True})
    if settings.HARVEST_ACCESS_TOKEN is None:
        return JsonResponse({"error": HARVEST_NOT_SET_UP}, status=400)

    try:
        harvest_project = client.get(f"projects/{project.harvest_id}")
    except requests.HTTPError as error:
        if error.response.status_code == 404:
            message = f"Project {project.harvest_id} is no longer in Harvest."
            return JsonResponse({"error": message}, status=404)
        return harvest_unavailable(error)
    except requests.RequestException as error:
        return harvest_unavailable(error)

    project.name = harvest_project["name"]
    project.harvest_notes = harvest_project["notes"]
    project.save()
    return JsonResponse({"success": True})
