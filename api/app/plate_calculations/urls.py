from django.urls import path

from plate_calculations.views import (
    correct_plate_background,
    delete_plate_calculation,
    log10_plate_measurement,
    normalize_plate_measurement,
)

urlpatterns = [
    path(
        "plates/<int:plate_id>/background_correction/",
        correct_plate_background,
        name="correct_plate_background",
    ),
    path(
        "plates/<int:plate_id>/log10/",
        log10_plate_measurement,
        name="log10_plate_measurement",
    ),
    path(
        "plates/<int:plate_id>/normalization/",
        normalize_plate_measurement,
        name="normalize_plate_measurement",
    ),
    path(
        "plates/<int:plate_id>/delete/",
        delete_plate_calculation,
        name="delete_plate_calculation",
    ),
]
