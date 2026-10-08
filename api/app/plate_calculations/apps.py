from django.apps import AppConfig


class PlateCalculationsConfig(AppConfig):
    """
    Calculations with the measurements of one plate, started from the plate page:
    background correction, log10 and %Activity. Each result is saved as a new
    measurement of the plate.
    """

    default_auto_field = "django.db.models.BigAutoField"
    name = "plate_calculations"
