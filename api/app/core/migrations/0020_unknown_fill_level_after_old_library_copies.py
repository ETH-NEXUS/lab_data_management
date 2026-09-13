"""
Turns the fill levels written by the old library plate copy into "unknown".

Library plates are imported without a volume. Until 2026-09 the copy action nevertheless
wrote a remaining fill level of 0 uL and a DMSO share of 100 % into every source well
(withdrawing 0 nL). Those zeros were never measured, but since zero now means "empty",
Recalculate would mark all those wells as problematic. The real fill level arrives with
the first Echo import, so the honest value until then is None ("not reported").

Only rows with all signs of the old copy are changed: 0 uL, 100 % DMSO, 0 nL withdrawn,
and a target well that was either deleted or lies on a copied plate ("_COPY").
Real Echo data, such as failed transfers (0 uL / 0 %), is left alone.
There is no way back: after the change these rows cannot be told apart from other
unreported values, and the zeros were wrong anyway.
"""

from django.db import migrations
from django.db.models import Q


def mark_old_copy_zeros_as_unknown(apps, schema_editor):
    WellWithdrawal = apps.get_model("core", "WellWithdrawal")
    old_copy_rows = WellWithdrawal.objects.filter(
        current_amount=0, current_dmso=100, amount=0
    ).filter(
        Q(target_well__isnull=True) | Q(target_well__plate__barcode__contains="_COPY")
    )
    updated = old_copy_rows.update(current_amount=None, current_dmso=None)
    print(f"\n  Marked {updated} fill levels from old library plate copies as unknown.")


class Migration(migrations.Migration):
    dependencies = [
        ("core", "0019_alter_threshold_amount_alter_threshold_dmso"),
    ]

    operations = [
        migrations.RunPython(mark_old_copy_zeros_as_unknown, migrations.RunPython.noop),
    ]
