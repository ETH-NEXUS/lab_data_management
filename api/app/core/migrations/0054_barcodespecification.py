# Generated by Django 4.1.7 on 2023-02-16 17:48

import core.basemodels
import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0053_measurementmetadata_measurement_identifier_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='BarcodeSpecification',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('prefix', models.CharField(max_length=100)),
                ('number_of_plates', models.IntegerField()),
                ('sides', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=20), size=None)),
                ('experiment', models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='barcode_specifications', to='core.experiment')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]