from django.urls import path

from background_correction.views import (
    correct_plate_background,
    log10_plate_measurement,
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
]
