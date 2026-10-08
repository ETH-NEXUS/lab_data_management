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
    Refused while measurements calculated from it are on the plate, e.g.
    "Lum1_log10_bc_R_median" when "Lum1_log10" is deleted: without it the plate
    page could no longer tell how they were calculated, nor offer to delete them.
    """
    with transaction.atomic():
        # The same lock as a calculation, so a calculation of this measurement
        # running at the same time is not mixed with its deletion
        Plate.objects.select_for_update(no_key=True).get(id=plate.id)
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
        # A calculation names its result after its source plus a suffix
        # (_log10, _bc_..., _inhibition_..., _activity_...)
        calculated_from_it = sorted(
            set(
                Measurement.objects.filter(
                    well__plate=plate, label__startswith=f"{label}_"
                ).values_list("label", flat=True)
            )
        )
        if calculated_from_it:
            raise ValidationError(
                f'{", ".join(calculated_from_it)} of the plate {plate.barcode} '
                f'were calculated from "{label}": delete them first.'
            )
        deleted, _ = measurements.delete()

    refresh_plate_views()
    return deleted
