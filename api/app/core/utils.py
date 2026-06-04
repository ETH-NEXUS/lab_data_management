from collections.abc import Iterable
from django.conf import settings
from django.db import transaction
from .models import Plate, PlateDetail, Well, WellCompound, WellDetail, WellWithdrawal


def build_copied_plate_barcode(source_barcode: str) -> str:
    """
    Build a unique barcode for a copied plate.
    Example input:
    {"barcode": "LIB_001"}
    Example output:
    {"barcode": "LIB_001_COPY_29.05.26"}
    """
    copy_suffix = f"_COPY"
    barcode = f"{source_barcode}{copy_suffix}"
    max_length = Plate._meta.get_field("barcode").max_length

    while True:
        if len(barcode) > max_length:
            raise ValueError(
                f"Copied barcode is longer than the allowed {max_length} characters."
            )
        if not Plate.objects.filter(barcode=barcode).exists():
            return barcode
        barcode = f"{barcode}{copy_suffix}"


def build_copy_withdrawal_metadata(source_well: Well, target_volume: float) -> dict:
    """
    Build withdrawal metadata for library plate copies.
    Example input:
    {"source_amount": 100, "target_volume": 20}
    Example output:
    {"current_amount": 80, "current_dmso": 100}
    """
    source_current_info = source_well.current_info
    remaining_amount = round(
        source_well.amount - target_volume,
        settings.FLOAT_PRECISION,
    )

    if source_current_info and source_current_info["current_dmso"] is not None:
        current_dmso = source_current_info["current_dmso"]
    else:
        current_dmso = 100

    return {
        "current_amount": max(remaining_amount, 0),
        "current_dmso": current_dmso,
    }


def copy_library_plates(
    source_plates: Iterable[Plate],
    target_volume: float,
) -> list[Plate]:
    """
    Copy library plates and create donor links back to the source wells.
    Copied data example:
    {
        "plate": {"barcode": "LIB_001", "library_id": 7},
        "wells": [{"position": 0, "status": "filled"}],
    }

    Returned data example:
    [{"barcode": "LIB_001_COPY_29.05.26", "library_id": 7}]
    """
    if target_volume <= 0:
        raise ValueError("Target volume must be greater than 0.")

    plates = list(source_plates)
    for source_plate in plates:
        if source_plate.library_id is None or source_plate.is_control_plate:
            raise ValueError("Only library plate copying is possible.")

    copied_plates = []
    with transaction.atomic():
        for source_plate in plates:
            copied_plate = Plate.objects.create(
                barcode=build_copied_plate_barcode(source_plate.barcode),
                dimension=source_plate.dimension,
                library=source_plate.library,
                archived=source_plate.archived,
                status=source_plate.status,
                use_as_template_to_select=source_plate.use_as_template_to_select,
            )

            source_wells = (
                source_plate.wells.select_related("sample", "type")
                .prefetch_related("well_compounds__compound")
                .order_by("position")
            )
            for source_well in source_wells:
                copied_well = Well.objects.create(
                    plate=copied_plate,
                    position=source_well.position,
                    sample=source_well.sample,
                    type=source_well.type,
                    status=source_well.status,
                    is_invalid=source_well.is_invalid,
                )

                source_well_compounds = list(source_well.well_compounds.all())
                if not source_well_compounds:
                    continue

                if source_well.amount < target_volume:
                    raise ValueError(
                        f"Well {source_well.hr_position} on plate {source_plate.barcode} "
                        f"does not have enough volume for {target_volume} nL."
                    )

                source_total_amount = source_well.initial_amount
                for source_well_compound in source_well_compounds:
                    if source_total_amount > 0:
                        copied_amount = round(
                            target_volume
                            * source_well_compound.amount
                            / source_total_amount,
                            settings.FLOAT_PRECISION,
                        )
                    else:
                        copied_amount = 0

                    WellCompound.objects.create(
                        well=copied_well,
                        compound=source_well_compound.compound,
                        amount=copied_amount,
                    )

                withdrawal_metadata = build_copy_withdrawal_metadata(
                    source_well,
                    target_volume,
                )
                WellWithdrawal.objects.create(
                    well=source_well,
                    target_well=copied_well,
                    amount=target_volume,
                    current_amount=withdrawal_metadata["current_amount"],
                    current_dmso=withdrawal_metadata["current_dmso"],
                )

            copied_plates.append(copied_plate)

    if copied_plates:
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    return copied_plates
