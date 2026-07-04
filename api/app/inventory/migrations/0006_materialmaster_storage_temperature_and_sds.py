from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("inventory", "0005_inventorystock_additional_sectors"),
    ]

    operations = [
        migrations.AddField(
            model_name="materialmaster",
            name="safety_data_sheet",
            field=models.FileField(
                blank=True,
                help_text="Optional safety data sheet attachment for this material.",
                null=True,
                upload_to="inventory/safety_data_sheets/",
            ),
        ),
        migrations.AddField(
            model_name="materialmaster",
            name="storage_temperature",
            field=models.CharField(
                blank=True,
                choices=[
                    ("4°C", "4°C"),
                    ("RT", "RT"),
                    ("-20°C", "-20°C"),
                    ("-80°C", "-80°C"),
                    ("LN", "LN"),
                ],
                help_text="Optional reagent storage temperature, e.g. 4°C, RT, -20°C, -80°C, or LN.",
                max_length=10,
                null=True,
            ),
        ),
    ]
