"""
Withdrawals that were counted twice.

Deleting the plates of an experiment and mapping the same Echo report again is a
common way to repair a mapping. The withdrawals of the library wells are kept
when their target plate is deleted (`target_well` becomes empty), so the second
mapping adds the same withdrawals again.

A withdrawal without target well can also be a real one: the liquid was
transferred, and the experiment plate was deleted later for another reason.
Two different Echo runs can also agree by chance in a single well. So a
mapping only counts as repeated when it repeats a whole plate
(see `repeated_mappings`).
"""

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from core.models import Plate, WellWithdrawal


@dataclass
class RepeatedMapping:
    """
    A plate mapped anew, and the withdrawals without target well it repeats, e.g.
    RepeatedMapping(source_plate=<Plate LIB_001>, target_plate=<Plate EXP_1>,
                    duplicates=[<WellWithdrawal A1 (20.0)>, ...])
    """

    source_plate: Plate
    target_plate: Plate
    duplicates: list = field(default_factory=list)


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


def repeated_mappings() -> list[RepeatedMapping]:
    """
    The plates that repeat a mapping whose target plate was deleted.

    A target plate repeats earlier withdrawals without target well when
    - it came from a report (a PlateMapping from the source plate; a plate copy
      has none) and was created after these withdrawals,
    - every withdrawal from the source plate into it has an Echo reading, and
    - each of them has a withdrawal without target well with the same well,
      amount and reading. One report gives the same numbers again; two
      different runs can agree in a well, but not in a whole plate.

    If the same report was mapped several times, every full set is a duplicate.
    A withdrawal is used for one repeated mapping only.
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
        source_plate = pool[0].well.plate
        target_plates = (
            Plate.objects.filter(
                mapped_from_plates__source_plate_id=source_plate_id,
                created_at__gt=pool[0].created_at,
            )
            .distinct()
            .order_by("created_at", "id")
        )
        for target_plate in target_plates:
            duplicates = repeated_withdrawals(pool, source_plate, target_plate)
            if not duplicates:
                continue
            for duplicate in duplicates:
                pool.remove(duplicate)
            results.append(RepeatedMapping(source_plate, target_plate, duplicates))
    return results


def repeated_withdrawals(pool: list, source_plate: Plate, target_plate: Plate) -> list:
    """
    The withdrawals of `pool` (without target well, of the source plate) that
    `target_plate` repeats as a whole, or [] if it does not.
    """
    mapped = list(
        WellWithdrawal.objects.filter(
            well__plate=source_plate, target_well__plate=target_plate
        )
    )
    if not mapped or not all(has_reading(withdrawal) for withdrawal in mapped):
        return []
    needed = Counter(withdrawal_key(withdrawal) for withdrawal in mapped)

    earlier_by_key = defaultdict(list)
    for orphan in pool:
        if orphan.created_at < target_plate.created_at:
            earlier_by_key[withdrawal_key(orphan)].append(orphan)

    # How many times the whole plate is repeated; 0 if one withdrawal is missing
    full_sets = min(len(earlier_by_key[key]) // count for key, count in needed.items())

    duplicates = []
    for key, count in needed.items():
        duplicates.extend(earlier_by_key[key][: full_sets * count])
    return duplicates
