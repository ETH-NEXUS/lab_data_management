"""
The Celery app of LDM. It runs tasks in the background, in the `celery` container.

A task is a function with @shared_task in a `tasks.py` of a Django app; the
worker finds these files by itself. Settings starting with CELERY_ in
ldm/settings.py configure the app (e.g. CELERY_BROKER_URL).
"""

import os

from celery import Celery

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ldm.settings")

app = Celery("ldm")
app.config_from_object("django.conf:settings", namespace="CELERY")
app.autodiscover_tasks()
