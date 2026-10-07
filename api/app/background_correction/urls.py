from django.urls import path

from background_correction.views import (
    correct_plate_background,
    log10_plate_measurement,
    percent_activity_of_plate,
)

urlpatterns = [
    path(
        "plates/<int:plate_id>/",
        correct_plate_background,
        name="correct_plate_background",
    ),
    path(
        "plates/<int:plate_id>/log10/",
        log10_plate_measurement,
        name="log10_plate_measurement",
    ),
    path(
        "plates/<int:plate_id>/activity/",
        percent_activity_of_plate,
        name="percent_activity_of_plate",
    ),
]
