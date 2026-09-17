"""
Small value conversions used by the mappers.
"""

import re
from datetime import datetime

from helpers.logger import logger

# The date and time formats of C10 files, found by their shape.
# Example: "10/14/2024 12:45:28" in a .txt file, "241014 125455" in a file name.
C10_DATETIME_FORMATS = {
    r"\d{1,2}/\d{1,2}/\d{4} \d{1,2}:\d{2}:\d{2}": "%m/%d/%Y %H:%M:%S",
    r"\d{6} \d{6}": "%y%m%d %H%M%S",
    r"\d{8} \d{6}": "%Y%m%d %H%M%S",
}


def parse_c10_datetime(date: str, time: str) -> datetime | None:
    """
    ("10/14/2024", "12:45:28") -> datetime(2024, 10, 14, 12, 45, 28).

    The result has no time zone: it is the local time of the instrument, like
    the other measurement dates. None if the date or time has an unknown format.
    """
    text = f"{date} {time}"
    for shape, date_format in C10_DATETIME_FORMATS.items():
        if re.fullmatch(shape, text):
            try:
                return datetime.strptime(text, date_format)
            except ValueError:
                # The right shape, but not a real date, e.g. "13/45/2024"
                return None
    return None


def convert_sci_to_float(sci_str):
    try:
        return float(sci_str)
    except ValueError:
        logger.error(f"Cannot convert {sci_str} to float")
        return None
