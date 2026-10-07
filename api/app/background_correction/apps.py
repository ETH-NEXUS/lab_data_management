from django.apps import AppConfig


class BackgroundCorrectionConfig(AppConfig):
    """
    The background correction of a plate: the median or mean of its reference
    wells is subtracted from the other wells, saved as a new measurement.
    """

    default_auto_field = "django.db.models.BigAutoField"
    name = "background_correction"
