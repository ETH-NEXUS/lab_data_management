from django.urls import path

from analysis.views import (
    download_analysis_result,
    list_analysis_results,
    start_analysis,
)

urlpatterns = [
    path("start/", start_analysis, name="start_analysis"),
    path("results/", list_analysis_results, name="list_analysis_results"),
    path("download/", download_analysis_result, name="download_analysis_result"),
]
