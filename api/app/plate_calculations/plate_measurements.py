"""
Reading and saving the measurements of one plate, shared by the calculations of
this app (background correction, log10, %Activity).
"""

from django.db import transaction
from rest_framework.exceptions import ValidationError

from core.models import (
    ExperimentDetail,
    Measurement,
    Plate,
    PlateDetail,
    WellDetail,
)

LABEL_MAX_LENGTH = Measurement._meta.get_field("label").max_length


def check_label_length(new_label: str) -> None:
    if len(new_label) > LABEL_MAX_LENGTH:
        raise ValidationError(
            f'The name of the new measurement "{new_label}" is longer than '
            f"{LABEL_MAX_LENGTH} characters."
        )


def values_by_time_point(plate: Plate, label: str) -> dict:
    """
    The values of one measurement of the plate, per time point and well id.
    A plate without the measurement is refused.

    {datetime(2025, 5, 16, 10, 0): {11: 10.0, 12: 4.0}}
    """
    measurements = Measurement.objects.filter(well__plate=plate, label=label)
    values: dict = {}
    for measurement in measurements:
        if measurement.measured_at not in values:
            values[measurement.measured_at] = {}
        values[measurement.measured_at][measurement.well_id] = measurement.value
    if not values:
        raise ValidationError(
            f'The plate {plate.barcode} has no measurement "{label}".'
        )
    return values


def well_ids_of_type(plate: Plate, well_type: str) -> set[int]:
    """The ids of the wells of one type on the plate, e.g. of the "R" wells."""
    return set(plate.wells.filter(type__name=well_type).values_list("id", flat=True))


def values_of_wells(
    plate: Plate,
    label: str,
    well_type: str,
    well_ids: set[int],
    well_values: dict[int, float],
    measured_at,
) -> list[float]:
    """
    The values of one time point in the wells of one type. A calculation that
    needs them (e.g. their median) is refused if there are none.

    well_values = {11: 10.0, 12: 4.0, 13: 7.0}, well_ids = {11, 13} -> [10.0, 7.0]
    """
    values = [well_values[well_id] for well_id in well_ids if well_id in well_values]
    if not values:
        raise ValidationError(
            f'The plate {plate.barcode} has no "{well_type}" wells with a '
            f'value of "{label}" measured at {measured_at}.'
        )
    return values


def replace_measurements(
    plate: Plate, new_label: str, new_measurements: list[Measurement]
) -> None:
    """
    Saves the new measurement of the plate instead of an earlier calculation of the
    same label, and refreshes the views the plate page reads. Refused: a result
    without a single value (it would only delete the earlier one), and a label of
    a measurement imported from a file (its values are not calculated here).
    """
    if not new_measurements:
        raise ValidationError(
            f'No well of the plate {plate.barcode} gets a value of "{new_label}".'
        )
    earlier = Measurement.objects.filter(well__plate=plate, label=new_label)
    # Only an import links its values to a file
    if earlier.filter(measurement_assignment__isnull=False).exists():
        raise ValidationError(
            f'The plate {plate.barcode} has an imported measurement "{new_label}", '
            "which a calculation does not replace."
        )

    with transaction.atomic():
        earlier.delete()
        Measurement.objects.bulk_create(new_measurements)

    PlateDetail.refresh(concurrently=True)
    WellDetail.refresh(concurrently=True)
    ExperimentDetail.refresh(concurrently=True)
