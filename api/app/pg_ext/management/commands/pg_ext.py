from django.core.management.base import BaseCommand
from django.core.management import call_command
import traceback
from colorful_logger import logger as log


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("action", type=str, help="The action to execute: apply")

    def apply(self):
        call_command("migrate", "--fake", "pg_ext", "zero")
        call_command("migrate")

    def handle(self, *args, **options):
        try:
            if options.get("action") == "apply":
                self.apply()
        except Exception as ex:
            log.error(ex)
            traceback.print_exc()
