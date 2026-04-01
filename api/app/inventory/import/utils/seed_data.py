import json
from pathlib import Path

from helpers.logger import logger
from .parsing import clean_value


def load_json_file(file_path):
    """
    Load JSON from file path.
    """
    logger.debug(f"Loading JSON seed file: {file_path}")
    with Path(file_path).open("r", encoding="utf-8") as file:
        return json.load(file)


def bulk_seed_named_values(model_class, values):
    """
    Create lookup values if they do not exist.
    Returns number of newly created records.
    """
    created_count = 0

    for value in values:
        value = clean_value(value)
        if not value:
            continue

        _, created = model_class.objects.get_or_create(name=value)
        if created:
            created_count += 1

    logger.debug(f"Seeded {created_count} new {model_class.__name__} records.")
    return created_count
