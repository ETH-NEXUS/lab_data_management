"""
Paths that come from the browser are checked here.

A page sends the path of a file it wants to read, write, delete or import.
Only paths inside a known folder are allowed, so a request cannot reach the
rest of the server.
"""

import os

from rest_framework.exceptions import ValidationError


def path_inside(path: str, root: str) -> str:
    """
    The path, if it is inside `root`, e.g. path_inside("/data/a.csv", "/data").

    ".." is resolved first, so "/data/../etc/passwd" is refused with an error
    that the page shows. A symbolic link inside the folder is allowed: the lab
    shares are mounted and linked in ways this code does not need to know.
    """
    if not isinstance(path, str) or not path:
        raise ValidationError("No path was given.")

    root = os.path.abspath(root)
    full_path = os.path.abspath(path)
    if full_path != root and not full_path.startswith(root + os.sep):
        raise ValidationError(f"The path must be inside {root}: {path}")
    return full_path
