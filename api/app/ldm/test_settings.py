# flake8: noqa
from .test_guard import ensure_tests_are_allowed

# Checked before anything else, so no test database is ever created outside development.
ensure_tests_are_allowed()

from .settings import *

# The tests run against the same Postgres server as the application, in a
# separate database that Django creates and drops per run. Sqlite cannot be
# used here: the migrations and the materialized views are Postgres specific.
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"

# The tests do not need the Redis container: they keep the cache in memory.
CACHES = {"default": {"BACKEND": "django.core.cache.backends.locmem.LocMemCache"}}

# Celery tasks run right away in the test process instead of in the celery container
CELERY_TASK_ALWAYS_EAGER = True
