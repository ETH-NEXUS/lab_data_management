# Generated by Django 4.1.7 on 2023-03-29 10:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0023_alter_measurement_unique_together_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='measurementfeature',
            name='abbrev',
            field=models.CharField(max_length=255, unique=True),
        ),
    ]