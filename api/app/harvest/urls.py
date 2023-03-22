from django.urls import path
from . import views

urlpatterns = [
    path("harvest_projects/", views.harvest_projects, name="harvest_projects"),
    path(
        "update_harvest_info/<int:project_id>/",
        views.update_harvest_info,
        name="update_harvest_info",
    ),
]
