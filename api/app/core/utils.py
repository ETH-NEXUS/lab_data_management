from collections.abc import Iterable
from django.db import transaction
from django.utils import timezone
from .models import Plate, PlateDetail, Well, WellCompound, WellDetail


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


def copy_library_plates(source_plates: Iterable[Plate]) -> list[Plate]:
    """
    Copy library plates without measurements or transfer history.
    Copied data example:
    {
        "plate": {"barcode": "LIB_001", "library_id": 7},
        "wells": [{"position": 0, "status": "filled"}],
    }

    Returned data example:
    [{"barcode": "LIB_001_COPY_29.05.26", "library_id": 7}]
    """
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

            for source_well in source_plate.wells.select_related(
                "sample", "type"
            ).prefetch_related("well_compounds__compound").order_by("position"):
                copied_well = Well.objects.create(
                    plate=copied_plate,
                    position=source_well.position,
                    sample=source_well.sample,
                    type=source_well.type,
                    status=source_well.status,
                    is_invalid=source_well.is_invalid,
                )
                for source_well_compound in source_well.well_compounds.all():
                    WellCompound.objects.create(
                        well=copied_well,
                        compound=source_well_compound.compound,
                        amount=source_well_compound.amount,
                    )

            copied_plates.append(copied_plate)

    if copied_plates:
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    return copied_plates
