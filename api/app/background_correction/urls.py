from django.urls import path

from background_correction.views import correct_plate_background

urlpatterns = [
    path(
        "plates/<int:plate_id>/",
        correct_plate_background,
        name="correct_plate_background",
    ),
]
