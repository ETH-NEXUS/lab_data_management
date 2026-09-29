from django.apps import AppConfig


class AnalysisConfig(AppConfig):
    """The statistical analysis of an experiment (the R reports in api/app/statistics)."""

    default_auto_field = "django.db.models.BigAutoField"
    name = "analysis"

    def ready(self):
        # Connects the signal that deletes the results of a deleted experiment
        from analysis import signals  # noqa: F401
