# Generated by Django 4.2.3 on 2023-09-12 14:54

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0007_well_is_invalid'),
    ]

    operations = [
        migrations.AlterField(
            model_name='measurement',
            name='label',
            field=models.CharField(default='none', max_length=50),
        ),
    ]
