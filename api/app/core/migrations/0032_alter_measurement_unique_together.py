# Generated by Django 4.1.4 on 2022-12-23 10:11

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0031_measurementfeature_remove_measurement_abbrev_and_more'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='measurement',
            unique_together={('well', 'feature')},
        ),
    ]
