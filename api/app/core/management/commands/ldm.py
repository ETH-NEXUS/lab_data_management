from django.core.management.base import BaseCommand
import logging

from core.models import PlateDetail, WellDetail, ExperimentDetail

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("action", type=str, help="The action to execute")

    def refresh_mat_views(self):
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

    def handle(self, *args, **options):
        try:
            if options.get("action") == "refresh_mat_views":
                self.refresh_mat_views()
        except Exception:
            logger.exception("ldm %s failed", options.get("action"))
