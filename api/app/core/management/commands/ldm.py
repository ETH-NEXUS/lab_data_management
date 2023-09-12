from django.core.management.base import BaseCommand
import traceback

from core.models import PlateDetail, WellDetail, ExperimentDetail


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
        except Exception as ex:
            print(ex)
            traceback.print_exc()
