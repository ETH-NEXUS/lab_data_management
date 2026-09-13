"""
Turns the fill levels written by library plate copies that were renamed later into "unknown".

Migration 0020 recognises a copy by a deleted target or a "_COPY" barcode. The lab renames
most copies to their real barcodes (for example Drug01_E -> Drug01_K), so their rows kept
the 0 uL that the old copy action wrote without any measurement.

Here a copy is recognised by how it is built instead of by its name: nothing was withdrawn
(0 nL), the fill level is 0 uL, the target well lies on a plate of the same library at the
same position, and there is no plate mapping record between the two plates. Echo imports
and csv mappings create such a record, so data from the instrument is left alone.
As in 0020 there is no way back, and the zeros were never measured.
"""

from django.db import migrations
from django.db.models import Exists, F, OuterRef


def mark_renamed_copy_zeros_as_unknown(apps, schema_editor):
    WellWithdrawal = apps.get_model("core", "WellWithdrawal")
    PlateMapping = apps.get_model("core", "PlateMapping")

    mapping_between_the_plates = PlateMapping.objects.filter(
        source_plate_id=OuterRef("well__plate_id"),
        target_plate_id=OuterRef("target_well__plate_id"),
    )
    renamed_copy_rows = WellWithdrawal.objects.filter(
        current_amount=0,
        amount=0,
        well__plate__library__isnull=False,
        target_well__isnull=False,
        target_well__plate__library_id=F("well__plate__library_id"),
        target_well__position=F("well__position"),
    ).filter(~Exists(mapping_between_the_plates))

    updated = renamed_copy_rows.update(current_amount=None, current_dmso=None)
    print(
        f"\n  Marked {updated} fill levels from renamed library plate copies as unknown."
    )


class Migration(migrations.Migration):
    dependencies = [
        ("core", "0020_unknown_fill_level_after_old_library_copies"),
    ]

    operations = [
        migrations.RunPython(
            mark_renamed_copy_zeros_as_unknown, migrations.RunPython.noop
        ),
    ]
