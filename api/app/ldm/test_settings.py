# flake8: noqa
from .test_guard import ensure_tests_are_allowed

# Checked before anything else, so no test database is ever created outside development.
ensure_tests_are_allowed()

from .settings import *

# The tests run against the same Postgres server as the application, in a
# separate database that Django creates and drops per run. Sqlite cannot be
# used here: the migrations and the materialized views are Postgres specific.
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"
