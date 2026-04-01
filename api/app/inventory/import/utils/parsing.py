from decimal import Decimal, InvalidOperation


def clean_value(value):
    """
    Normalize empty / whitespace-only values to None.
    """
    if value is None:
        return None

    value = str(value).strip()

    if value == "":
        return None

    if value.lower() in {"nan", "none", "null"}:
        return None

    return value


def to_decimal(value, default=None):
    """
    Parse a decimal safely.

    Examples:
    - 10
    - 0.8
    - 1.366666667
    """
    value = clean_value(value)
    if value is None:
        return default

    try:
        return Decimal(value)
    except (InvalidOperation, ValueError):
        return default


def to_bool(value):
    """
    Convert common truthy values to bool.
    """
    value = clean_value(value)
    if not value:
        return False

    return value.lower() in {"1", "true", "yes", "y", "x", "favorite", "favourite"}


def merge_notes(old_notes, new_notes):
    """
    Merge notes without losing existing text.
    """
    old_notes = clean_value(old_notes)
    new_notes = clean_value(new_notes)

    if old_notes and new_notes:
        if new_notes in old_notes:
            return old_notes
        return f"{old_notes}\n{new_notes}"

    return old_notes or new_notes


def split_attributes(item_property_text):
    """
    Convert Excel property strings into reusable attribute names.

    Examples:
    - 'sterile & filtered' -> ['sterile', 'filtered']
    - 'pure' -> ['pure']
    """
    if not item_property_text:
        return []

    normalized = item_property_text.replace("&", ",")
    parts = [part.strip() for part in normalized.split(",")]

    cleaned_parts = []
    for part in parts:
        if not part:
            continue
        cleaned_parts.append(part)

    return cleaned_parts