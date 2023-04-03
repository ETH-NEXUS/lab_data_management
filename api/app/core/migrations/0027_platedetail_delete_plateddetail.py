# Generated by Django 4.1.7 on 2023-03-31 14:23

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("core", "0026_plateddetail_welldetail_and_more"),
    ]

    operations = [
        migrations.CreateModel(
            name="PlateDetail",
            fields=[
                ("id", models.BigIntegerField(primary_key=True, serialize=False)),
                ("barcode", models.CharField(max_length=50)),
                ("dimension", models.JSONField()),
                (
                    "measurement_labels",
                    django.contrib.postgres.fields.ArrayField(
                        base_field=models.TextField(blank=True, null=True), size=None
                    ),
                ),
                ("min_max", models.JSONField()),
                ("stats", models.JSONField()),
            ],
            options={
                "db_table": "core_platedetail",
                "managed": False,
            },
        ),
        migrations.DeleteModel(
            name="PlatedDetail",
        ),
    ]
