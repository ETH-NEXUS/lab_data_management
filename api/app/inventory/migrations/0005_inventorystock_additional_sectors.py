from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("inventory", "0004_inventorystocktablepreference"),
    ]

    operations = [
        migrations.AddField(
            model_name="inventorystock",
            name="additional_sectors",
            field=models.ManyToManyField(
                blank=True,
                help_text=(
                    "Additional physical storage sectors for the same stock item. "
                    "Example: one consumable split across multiple shelves in the same room."
                ),
                related_name="secondary_stock_entries",
                to="inventory.sector",
            ),
        ),
    ]
