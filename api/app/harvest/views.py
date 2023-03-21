import json

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

from .harvest_client import HarvestClient
from django.conf import settings
from core.models import Project, HarvestProjectMapping
from django.shortcuts import get_object_or_404

client = HarvestClient(settings.HARVEST_ACCESS_TOKEN, settings.HARVEST_ACCOUNT_ID)


def harvest_projects(request):
    projects = client.get("projects")
    return JsonResponse(projects)


@csrf_exempt
def update_harvest_info(request, project_id):
    available_harvest_projects = client.get("projects")["projects"]
    project = get_object_or_404(Project, id=project_id)
    harvest_project_mapping = get_object_or_404(HarvestProjectMapping, project=project)
    harvest_project = list(
        filter(
            lambda x: x["id"] == harvest_project_mapping.harvest_id,
            available_harvest_projects,
        )
    )[0]

    if not harvest_project:
        return JsonResponse({"success": False, "error": "Harvest project not found."})

    harvest_project_mapping.name = project.name = harvest_project["name"]
    harvest_project_mapping.save()
    project.save()

    return JsonResponse({"success": True})
