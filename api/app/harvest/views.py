import json

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

from .harvest_client import HarvestClient
from django.conf import settings
from core.models import Project
from django.shortcuts import get_object_or_404
from os import environ

client = HarvestClient(settings.HARVEST_ACCESS_TOKEN, settings.HARVEST_ACCOUNT_ID)

HARVEST_PROJECT_FILTER = environ.get("HARVEST_PROJECT_FILTER")


def harvest_projects(request, filter_string=HARVEST_PROJECT_FILTER):
    if settings.HARVEST_ACCESS_TOKEN is None:
        return JsonResponse({"projects": []})
    projects = client.get("projects")
    if filter_string:
        projects = {
            "projects": list(
                filter(lambda x: filter_string in x["name"], projects["projects"])
            )
        }
    return JsonResponse(projects)


@csrf_exempt
def update_harvest_info(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    available_harvest_projects = client.get("projects")["projects"]

    if project.harvest_id:
        harvest_project = list(
            filter(lambda x: x["id"] == project.harvest_id, available_harvest_projects)
        )[0]
        project.name = harvest_project["name"]
        project.harvest_notes = harvest_project["notes"]
        project.save()

    return JsonResponse({"success": True})
