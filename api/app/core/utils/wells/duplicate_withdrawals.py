"""
Withdrawals that may have been counted twice.

Deleting the plates of an experiment and mapping the same Echo report again is a
common way to repair a mapping. The withdrawals of the library wells are kept
when their target plate is deleted (`target_well` becomes empty), so the second
mapping adds the same withdrawals again.

The database cannot tell this apart from a real transfer: the Echo writes one
report per destination plate, and replicate plates of the same library plate get
the same wells, volumes and readings (seen in 2024_snl_pruschy_radiation: four
plates from L3900-2_11 with identical numbers). If one of them is deleted, the
others repeat its withdrawals exactly. So this module only finds candidates;
someone who knows that a plate was mapped anew has to confirm it.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from core.models import Plate, WellWithdrawal


@dataclass
class RepeatedMapping:
    """
    A plate that repeats withdrawals without target well as a whole, e.g.
    RepeatedMapping(source_plate=<Plate LIB_001>, target_plate=<Plate EXP_1>,
                    one_set=[<WellWithdrawal A1 (20.0)>, ...], full_sets=1)

    `one_set` is what one mapping of the plate added; `full_sets` says how many
    such sets there are among the withdrawals without target well.
    """

    source_plate: Plate
    target_plate: Plate
    one_set: list = field(default_factory=list)
    full_sets: int = 0


def withdrawal_key(withdrawal: WellWithdrawal) -> tuple:
    """What the same transfer has in common, e.g. (well id, 20.0, 9.86, 96.1)."""
    return (
        withdrawal.well_id,
        withdrawal.amount,
        withdrawal.current_amount,
        withdrawal.current_dmso,
    )


def has_reading(withdrawal: WellWithdrawal) -> bool:
    """True if the Echo reported the volume and the DMSO of the source well."""
    return withdrawal.current_amount is not None and withdrawal.current_dmso is not None


def repeated_mappings(target_barcode: str | None = None) -> list[RepeatedMapping]:
    """
    The plates that repeat withdrawals without target well as a whole: each
    one may be a report mapped again, or a real transfer with the same numbers.

    A target plate counts when
    - it came from a report (a PlateMapping from the source plate; a plate copy
      has none) and was created after these withdrawals,
    - every withdrawal from the source plate into it has an Echo reading, and
    - each of them has a withdrawal without target well with the same well,
      amount and reading.

    With `target_barcode` only that plate is looked at.
    """
    orphans_by_source = defaultdict(list)
    orphans = (
        WellWithdrawal.objects.filter(target_well__isnull=True)
        .select_related("well__plate__dimension")
        .order_by("created_at", "id")
    )
    for orphan in orphans:
        if has_reading(orphan):
            orphans_by_source[orphan.well.plate_id].append(orphan)

    results = []
    for source_plate_id, pool in orphans_by_source.items():
        target_plates = Plate.objects.filter(
            mapped_from_plates__source_plate_id=source_plate_id,
            created_at__gt=pool[0].created_at,
        )
        if target_barcode is not None:
            target_plates = target_plates.filter(barcode=target_barcode)
        for target_plate in target_plates.distinct().order_by("created_at", "id"):
            repeated = repeated_mapping(pool, target_plate)
            if repeated is not None:
                results.append(repeated)
    return results


def repeated_mapping(pool: list, target_plate: Plate) -> RepeatedMapping | None:
    """
    How `target_plate` repeats the withdrawals of `pool` (without target well,
    all of one source plate), or None if it does not repeat them as a whole.
    """
    source_plate = pool[0].well.plate
    mapped = list(
        WellWithdrawal.objects.filter(
            well__plate=source_plate, target_well__plate=target_plate
        )
    )
    if not mapped or not all(has_reading(withdrawal) for withdrawal in mapped):
        return None
    needed = Counter(withdrawal_key(withdrawal) for withdrawal in mapped)

    earlier_by_key = defaultdict(list)
    for orphan in pool:
        if orphan.created_at < target_plate.created_at:
            earlier_by_key[withdrawal_key(orphan)].append(orphan)

    # How many times the whole plate is repeated; 0 if one withdrawal is missing
    full_sets = min(len(earlier_by_key[key]) // count for key, count in needed.items())
    if full_sets == 0:
        return None

    # One set, the newest withdrawals: a mapping that is repeated is most likely
    # the last one before the new plate. They are equal, so it does not change
    # any amount which of them is taken.
    one_set = []
    for key, count in needed.items():
        one_set.extend(earlier_by_key[key][-count:])
    return RepeatedMapping(source_plate, target_plate, one_set, full_sets)
