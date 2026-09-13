#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys

from ldm.test_guard import ensure_tests_are_allowed


def main():
    """Run administrative tasks."""
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        # The tests create a test database on the connected server, so they only
        # run in development and always with the test settings.
        ensure_tests_are_allowed()
        os.environ["DJANGO_SETTINGS_MODULE"] = "ldm.test_settings"
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ldm.settings")
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == "__main__":
    main()
