from collections.abc import Iterable
from django.conf import settings
from django.db import transaction
from core.models import (
    Plate,
    PlateDetail,
    Well,
    WellCompound,
    WellDetail,
    WellWithdrawal,
)
from core.utils.wells.volume_units import nanoliter_to_microliter


def build_copied_plate_barcode(source_barcode: str) -> str:
    """
    Build a unique barcode for a copied plate.
    Example input:
    {"barcode": "LIB_001"}
    Example output:
    {"barcode": "LIB_001_COPY_29.05.26"}
    """
    copy_suffix = "_COPY"
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


def build_copy_withdrawal_metadata(source_well: Well, withdrawal_volume: float) -> dict:
    """
    Build withdrawal metadata for library plate copies.
    The fill level of the source well is only known from what it reported in an
    earlier withdrawal, so the copy carries that level over, minus the copied
    volume. `withdrawal_volume` is given in nanoliter, `current_amount` is
    returned in microliter, the unit of the threshold on the messages page.
    Both values are None when the source well has never reported them, which
    means "unknown" and not "empty".
    Example input:
    {"source_well": {"current_info": {"current_amount": 9.5, "current_dmso": 95}},
     "withdrawal_volume": 500}
    Example output:
    {"current_amount": 9.0, "current_dmso": 95}
    """
    source_current_info = source_well.current_info or {}
    source_current_amount = source_current_info.get("current_amount")

    if source_current_amount is None:
        current_amount = None
    else:
        remaining_amount = source_current_amount - nanoliter_to_microliter(
            withdrawal_volume
        )
        current_amount = round(max(remaining_amount, 0), settings.FLOAT_PRECISION)

    return {
        "current_amount": current_amount,
        "current_dmso": source_current_info.get("current_dmso"),
    }


def get_copied_compound_amounts(
    source_well_compounds: list[WellCompound],
    source_total_amount: float,
    target_volume: float,
) -> list[tuple[WellCompound, float]]:
    """
    Calculate copied compound amounts for a library plate copy.
    Example input:
    {"source_total_amount": 0, "target_volume": 20, "compound_count": 1}
    Example output:
    [{"compound": "source compound", "amount": 20}]
    """
    if source_total_amount > 0:
        copied_amounts = []
        for source_well_compound in source_well_compounds:
            source_amount = source_well_compound.amount or 0
            copied_amount = round(
                target_volume * source_amount / source_total_amount,
                settings.FLOAT_PRECISION,
            )
            copied_amounts.append((source_well_compound, copied_amount))
        return copied_amounts

    copied_amounts = []
    remaining_volume = target_volume
    compound_count = len(source_well_compounds)
    for index, source_well_compound in enumerate(source_well_compounds):
        if index == compound_count - 1:
            copied_amount = round(remaining_volume, settings.FLOAT_PRECISION)
        else:
            copied_amount = round(
                target_volume / compound_count,
                settings.FLOAT_PRECISION,
            )
            remaining_volume -= copied_amount
        copied_amounts.append((source_well_compound, copied_amount))
    return copied_amounts


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
            # A copy is a fresh plate, so the automatic "running low" mark of the
            # source does not apply to it. A status set by hand is copied as before.
            copied_plate_status = source_plate.status
            if copied_plate_status == "empty_wells":
                copied_plate_status = None
            copied_plate = Plate.objects.create(
                barcode=build_copied_plate_barcode(source_plate.barcode),
                dimension=source_plate.dimension,
                library=source_plate.library,
                # A copy is a plate that exists in storage, so it is never archived,
                # even when the plate it was copied from is.
                archived=False,
                status=copied_plate_status,
                use_as_template_to_select=source_plate.use_as_template_to_select,
            )

            source_wells = (
                source_plate.wells.select_related("sample", "type")
                .prefetch_related("well_compounds__compound")
                .order_by("position")
            )
            for source_well in source_wells:
                # The same for the wells: a copied well is not empty.
                copied_well_status = source_well.status
                if copied_well_status == "empty":
                    copied_well_status = None
                copied_well = Well.objects.create(
                    plate=copied_plate,
                    position=source_well.position,
                    sample=source_well.sample,
                    type=source_well.type,
                    status=copied_well_status,
                    is_invalid=source_well.is_invalid,
                )

                source_well_compounds = list(source_well.well_compounds.all())
                if not source_well_compounds:
                    continue

                source_total_amount = source_well.initial_amount
                has_source_volume = source_total_amount > 0

                if has_source_volume and source_well.amount < target_volume:
                    raise ValueError(
                        f"Well {source_well.hr_position} on plate {source_plate.barcode} "
                        f"does not have enough volume for {target_volume} nL."
                    )

                copied_amounts = get_copied_compound_amounts(
                    source_well_compounds,
                    source_total_amount,
                    target_volume,
                )
                for source_well_compound, copied_amount in copied_amounts:
                    WellCompound.objects.create(
                        well=copied_well,
                        compound=source_well_compound.compound,
                        amount=copied_amount,
                    )

                # A library well without a volume (an SDF value like "<24", or a
                # library whose amounts were not filled with `fill_sdf_amounts`)
                # has no nanoliter bookkeeping to withdraw from. Withdrawing
                # anyway would push Well.amount below zero, which then shows up
                # as a negative volume on the plate page.
                withdrawal_volume = target_volume if has_source_volume else 0
                withdrawal_metadata = build_copy_withdrawal_metadata(
                    source_well,
                    withdrawal_volume,
                )
                WellWithdrawal.objects.create(
                    well=source_well,
                    target_well=copied_well,
                    amount=withdrawal_volume,
                    current_amount=withdrawal_metadata["current_amount"],
                    current_dmso=withdrawal_metadata["current_dmso"],
                )

            copied_plates.append(copied_plate)

    if copied_plates:
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    return copied_plates
