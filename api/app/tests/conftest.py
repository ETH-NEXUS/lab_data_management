"""
pytest loads this file before any test runs, whatever settings module is in use.

The guard is needed here and not only in ldm/test_settings.py: pytest takes
DJANGO_SETTINGS_MODULE from the environment first, and .env sets it to
ldm.settings, so the test settings, and the guard in them, are not loaded.
"""

from ldm.test_guard import ensure_tests_are_allowed

ensure_tests_are_allowed()
