# Generated by Django 4.1.7 on 2023-04-06 15:30

import core.models
import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0002_manual_plate_well_mat_view'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExperimentDetail',
            fields=[
                ('id', models.BigIntegerField(primary_key=True, serialize=False)),
                ('project_id', models.BigIntegerField(blank=True, null=True)),
                ('measurement_labels', django.contrib.postgres.fields.ArrayField(base_field=models.TextField(blank=True, null=True), size=None)),
                ('measurement_timestamps', core.models.DictField()),
                ('stats', core.models.DictField()),
                ('overall_stats', core.models.DictField()),
            ],
            options={
                'db_table': 'core_experimentdetail',
                'managed': False,
            },
        ),
    ]