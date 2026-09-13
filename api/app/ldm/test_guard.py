"""
Stops the tests from running anywhere but in the development environment.

The tests create and drop a test database on the database server the application
is connected to. On production that is the production server, so running them
there must not be possible, not even by accident.
Only the development containers set LDM_ALLOW_TESTS=1 (docker-compose.dev.yml).
"""

import os
import sys

ALLOW_TESTS_VARIABLE = "LDM_ALLOW_TESTS"


def ensure_tests_are_allowed():
    """
    Exit before Django connects to any database when tests are not allowed here.
    """
    if os.environ.get(ALLOW_TESTS_VARIABLE) == "1":
        return

    sys.exit(
        "Tests are disabled here: they would create a test database on this "
        f"database server. They only run where {ALLOW_TESTS_VARIABLE}=1 is set, "
        "which docker-compose.dev.yml does."
    )
