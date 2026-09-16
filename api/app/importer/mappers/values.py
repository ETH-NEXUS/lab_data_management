"""
Small value conversions used by the mappers.
"""

from datetime import datetime as dt
from datetime import timezone

from helpers.logger import logger

GLOBAL_NOW = dt.now(timezone.utc)


def convert_string_to_datetime(date_str, time_str):
    try:
        formatted_date_str = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
        formatted_time_str = f"{time_str[:2]}:{time_str[2:4]}:{time_str[4:]}"
        combined_str = f"{formatted_date_str} {formatted_time_str}"
        datetime_obj = dt.strptime(combined_str, "%Y-%m-%d %H:%M:%S")

        datetime_obj = datetime_obj.replace(tzinfo=timezone.utc)

        return datetime_obj.isoformat()
    except ValueError as e:

        logger.warning(
            f"Cannot convert {date_str} {time_str} to datetime: {e}. The current time will be used instead."
        )
        return GLOBAL_NOW.isoformat()


def convert_sci_to_float(sci_str):
    try:
        return float(sci_str)
    except ValueError:
        logger.error(f"Cannot convert {sci_str} to float")
        return None
