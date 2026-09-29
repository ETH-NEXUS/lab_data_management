"""
The rows of the "add experiment data" form: one row per plate and measurement
label, with the values that were saved before.
"""


def with_saved_values(new_rows: list[dict], saved_rows: list[dict]) -> list[dict]:
    """
    The rows of every plate and label, with the saved values where a row was saved.
    Before, only the saved rows were shown, so a label that was never saved could
    not be added any more. A saved row whose measurements are gone is kept.

    new_rows:   [{"plate_barcode": "SP_1", "measurement_label": "Lum", "condition": "", ...},
                 {"plate_barcode": "SP_1", "measurement_label": "Fluo", "condition": "", ...}]
    saved_rows: [{"plate_barcode": "SP_1", "measurement_label": "Lum", "condition": "KO", ...}]
    returns:    [{... "Lum", "condition": "KO", ...}, {... "Fluo", "condition": "", ...}]
    """
    saved_by_key = {
        (row["plate_barcode"], row["measurement_label"]): row for row in saved_rows
    }
    rows = []
    for row in new_rows:
        key = (row["plate_barcode"], row["measurement_label"])
        rows.append(saved_by_key.pop(key, row))
    rows.extend(saved_by_key.values())
    return rows
