"""
The two ways the M1000 names its files. The lab chooses the format on the
management page; the file name gives the plate barcode and, if the file itself
does not, the time of the measurement.

- "date_time_barcode" (the older setting): "20240610-121212_demo_1.asc",
  date and time may be left out: "demo_1.asc". The "Evaluated" files of the
  M1000 write the date with 6 digits, month first as in the newer setting:
  "093026-154654_RKS_300926_1.asc" (09/30/26)
- "barcode_date_time" (the newer setting): "RKS_300926_3_093026_165454.asc",
  barcode "RKS_300926_3", date 09/30/26 (month, day, year), time 16:54:54
"""

import os
import re
from datetime import datetime
from typing import TypedDict

from django.core.management.base import CommandError

from importer.mappers.base import SkipFile

DATE_TIME_BARCODE = "date_time_barcode"
BARCODE_DATE_TIME = "barcode_date_time"
FILE_NAME_FORMATS = (DATE_TIME_BARCODE, BARCODE_DATE_TIME)

FILE_NAME_PATTERNS = {
    DATE_TIME_BARCODE: r"^(?:(?P<date>[0-9]+)-(?P<time>[0-9]+)_)?(?P<barcode>[^\.]+)\.asc$",
    # Date and time are at the end, so the barcode may contain underscores ("RKS_300926_3")
    BARCODE_DATE_TIME: r"^(?P<barcode>.+)_(?P<date>[0-9]{6})_(?P<time>[0-9]{6})\.asc$",
}

# How the date and the time are written in the name, for datetime.strptime,
# by the number of digits of the date
DATE_TIME_FORMATS = {
    DATE_TIME_BARCODE: {
        8: "%Y%m%d %H%M%S",  # "20240610 121212"
        6: "%m%d%y %H%M%S",  # "093026 154654" of the Evaluated files, month first
    },
    BARCODE_DATE_TIME: {
        6: "%m%d%y %H%M%S",  # "093026 165454", the month comes first
    },
}

FILE_NAME_EXAMPLES = {
    DATE_TIME_BARCODE: "20240610-121212_demo_1.asc",
    BARCODE_DATE_TIME: "RKS_300926_3_093026_165454.asc",
}


class M1000FileName(TypedDict):
    barcode: str  # e.g. "demo_1"
    measured_at: datetime | None  # None if the name has no (readable) date


def read_file_name(path: str, file_name_format: str | None) -> M1000FileName:
    """
    ("/data/RKS_300926_3_093026_165454.asc", "barcode_date_time")
    -> {"barcode": "RKS_300926_3", "measured_at": datetime(2026, 9, 30, 16, 54, 54)}

    Without a format the older one is used, as before there was a choice.
    """
    chosen_format = file_name_format or DATE_TIME_BARCODE
    if chosen_format not in FILE_NAME_PATTERNS:
        raise CommandError(
            f"Unknown file name format '{chosen_format}'. "
            f"Known formats: {', '.join(FILE_NAME_FORMATS)}."
        )

    file_name = os.path.basename(path)
    match = re.match(FILE_NAME_PATTERNS[chosen_format], file_name)
    if not match:
        # Other files can lie in the folder, e.g. "30092026-001.asc" next to the
        # files of the plates; they are skipped, the others are still mapped
        raise SkipFile(
            f"its name does not match the chosen format {chosen_format}, "
            f"e.g. {FILE_NAME_EXAMPLES[chosen_format]}."
        )
    # The older pattern takes any name, so a newer name would become a wrong
    # barcode, e.g. "RKS_300926_3_093026_165454", and a new plate of that name
    if chosen_format == DATE_TIME_BARCODE and re.match(
        FILE_NAME_PATTERNS[BARCODE_DATE_TIME], file_name
    ):
        raise CommandError(
            f"The file name {file_name} looks like the format {BARCODE_DATE_TIME}, "
            f"e.g. {FILE_NAME_EXAMPLES[BARCODE_DATE_TIME]}. Choose that format."
        )

    return {
        "barcode": match.group("barcode"),
        "measured_at": date_from_match(
            match, DATE_TIME_FORMATS[chosen_format], file_name
        ),
    }


def date_from_match(
    match: re.Match, date_time_formats: dict[int, str], file_name: str
) -> datetime | None:
    """
    date "093026", time "165454", {6: "%m%d%y %H%M%S"} -> datetime(2026, 9, 30, 16, 54, 54)

    A date with another number of digits is refused: datetime.strptime would read
    "093026" with "%Y%m%d" as the year 930, month 2, day 6, without an error.
    """
    date = match.group("date")
    if not date:
        return None
    if len(date) not in date_time_formats:
        digits = " or ".join(str(length) for length in date_time_formats)
        raise CommandError(
            f"The date {date} in the file name {file_name} has {len(date)} digits, "
            f"but {digits} are expected (e.g. 20240610 or 093026, month first)."
        )
    try:
        return datetime.strptime(
            f"{date} {match.group('time')}", date_time_formats[len(date)]
        )
    except ValueError:
        return None
