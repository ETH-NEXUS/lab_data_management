# Generated by Django 4.1.7 on 2023-03-07 15:05

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0004_alter_measurementfeature_abbrev'),
    ]

    operations = [
        migrations.AlterField(
            model_name='measurementfeature',
            name='abbrev',
            field=models.CharField(max_length=20, unique=True),
        ),
    ]
