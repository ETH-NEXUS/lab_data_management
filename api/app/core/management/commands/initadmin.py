from django.core.management.base import BaseCommand
from django.contrib.auth import get_user_model
from os import environ
from helpers.logger import logger
from core.models import Threshold


class Command(BaseCommand):
    def handle(self, *args, **options):
        superuser = environ.get("DJANGO_SU_NAME", "admin")
        User = get_user_model()
        if not User.objects.filter(username=superuser).exists():
            User.objects.create_superuser(
                superuser,
                environ.get("DJANGO_SU_EMAIL", "admin@admin.com"),
                environ.get("DJANGO_SU_PASSWORD", "superuser"),
            )
            logger.debug("Created superuser account")
        else:
            logger.warning("Superuser exists")
        if not Threshold.objects.exists():
            Threshold.objects.create()
            logger.debug("Created threshold")
