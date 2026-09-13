"""
Application version, the bundled docs and refreshing the materialized views.
"""

import mimetypes
import os
import re
from os import environ
import traceback
from django.conf import settings
from django.http import Http404
from django.http import HttpResponse
from django.views.generic import View
from rest_framework import views
from django.http import JsonResponse
from rest_framework.response import Response
from ..models import (
    PlateDetail,
    WellDetail,
)
from ..serializers import ExperimentDetail


class VersionView(views.APIView):
    def get(self, request, format=None):
        return Response({"version": environ.get("GIT_VERSION", "N/A")})


class DocsView(View):
    def get(self, request, uri, **kwargs):
        docs_dir = os.path.join(settings.BASE_DIR, "docs", "site")
        if uri == "":
            uri = "index.html"
        file_path = os.path.join(docs_dir, uri)

        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "index.html")
        if not os.path.isfile(file_path):
            raise Http404("File not found")

        with open(file_path, "r", encoding="utf8", errors="replace") as f:
            content = f.read()

        # Replace all emojis with an empty string
        content = re.sub(r"[^\x00-\x7F]+", "", content)
        mime_type = mimetypes.guess_type(file_path)
        return HttpResponse(content, content_type=mime_type[0])


def refresh(request):
    """
    Hard refresh on the ui side.
    When we delete measurements in admin, it will not be imeediately reflected on the UI,
    because we need to refresh materialized views.
    This function is triggerd by the refresh button on the ui.
    """
    try:
        if request.method == "GET":
            PlateDetail.refresh(concurrently=True)
            WellDetail.refresh(concurrently=True)
            ExperimentDetail.refresh(concurrently=True)
            return JsonResponse({"status": "Data refreshed successfully"}, status=200)
    except Exception as e:
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)
