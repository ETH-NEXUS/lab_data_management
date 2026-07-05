from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("inventory", "0006_materialmaster_storage_temperature_and_sds"),
    ]

    operations = [
        migrations.AddField(
            model_name="inventorystock",
            name="source_order",
            field=models.ForeignKey(
                blank=True,
                help_text=(
                    "Optional order that this stock entry was created from. "
                    "Used to trace one stock item back to the originating order."
                ),
                null=True,
                on_delete=models.SET_NULL,
                related_name="created_stock_entries",
                to="inventory.order",
            ),
        ),
    ]
