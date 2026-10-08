"""
Deleting a measurement that was calculated on the plate page, e.g. "Lum1_log10".
A measurement imported from a file is never deleted here.
"""

from django.db import transaction
from rest_framework.exceptions import ValidationError

from core.models import Measurement, Plate
from plate_calculations.plate_measurements import refresh_plate_views


def delete_calculated_measurement(plate: Plate, label: str) -> int:
    """
    Deletes every value of one calculated measurement of the plate (all wells and
    reads) and returns how many values were deleted, e.g. 64.
    Measurements calculated from it stay (e.g. "Lum1_log10_bc_R_median" when
    "Lum1_log10" is deleted).
    """
    with transaction.atomic():
        # The same lock as a calculation, so a calculation of this measurement
        # running at the same time is not mixed with its deletion
        Plate.objects.select_for_update().get(id=plate.id)
        measurements = Measurement.objects.filter(well__plate=plate, label=label)
        if not measurements.exists():
            raise ValidationError(
                f'The plate {plate.barcode} has no measurement "{label}".'
            )
        # Only an import links its values to a file
        if measurements.filter(measurement_assignment__isnull=False).exists():
            raise ValidationError(
                f'"{label}" of the plate {plate.barcode} was imported from a file, '
                "so it cannot be deleted here."
            )
        deleted, _ = measurements.delete()

    refresh_plate_views()
    return deleted
