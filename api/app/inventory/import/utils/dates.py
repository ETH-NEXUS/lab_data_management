from datetime import datetime
from django.utils import timezone

from helpers.logger import logger
from .parsing import clean_value


def parse_flexible_date(value):
    """
    Parse many real-world date formats from the Excel export.

    Supports:
    - 2019-11
    - 2023-01-28
    - 27.07.2028
    - 20280610
    - 06/2030
    - Sep 28
    - Dez 28
    """
    value = clean_value(value)
    if not value:
        logger.debug("No date value provided; returning None.")
        return None

    replacements = {
        "Dez": "Dec",
        "Mär": "Mar",
        "Mai": "May",
        "Okt": "Oct",
    }

    for old, new in replacements.items():
        value = value.replace(old, new)

    date_formats = [
        "%Y-%m-%d",
        "%Y-%m",
        "%d.%m.%Y",
        "%d.%m.%y",
        "%d/%m/%Y",
        "%m/%Y",
        "%b %y",
        "%B %y",
        "%Y%m%d",
    ]

    for fmt in date_formats:
        try:
            parsed = datetime.strptime(value, fmt)
            if fmt in {"%Y-%m", "%m/%Y", "%b %y", "%B %y"}:
                return parsed.replace(day=1).date()
            return parsed.date()
        except ValueError:
            continue

    logger.debug(f"Could not parse date value: {value!r}")
    return None


def parse_flexible_datetime(value):
    """
    Parse a datetime. Fall back to 'now' if missing or invalid.
    """
    value = clean_value(value)
    if not value:
        logger.debug("No datetime value provided; using current time.")
        return timezone.now()

    datetime_formats = [
        "%Y-%m-%d",
        "%d.%m.%Y",
        "%d.%m.%y",
        "%Y-%m-%d %H:%M:%S",
    ]

    for fmt in datetime_formats:
        try:
            parsed = datetime.strptime(value, fmt)
            return timezone.make_aware(parsed, timezone.get_current_timezone())
        except ValueError:
            continue

    logger.debug(f"Could not parse datetime value: {value!r}; using current time.")
    return timezone.now()
