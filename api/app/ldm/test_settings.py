# flake8: noqa
from .settings import *

# The tests run against the same Postgres server as the application, in a
# separate database that Django creates and drops per run. Sqlite cannot be
# used here: the migrations and the materialized views are Postgres specific.
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"
