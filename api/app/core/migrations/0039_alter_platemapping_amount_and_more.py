# Generated by Django 4.1.4 on 2023-01-05 13:43

import django.core.validators
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0038_platemapping_amount_alter_measurement_well_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='platemapping',
            name='amount',
            field=models.FloatField(default=None, null=True, validators=[django.core.validators.MinValueValidator(0)]),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='amount_column',
            field=models.CharField(max_length=50, null=True),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='delimiter',
            field=models.CharField(default=',', max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='from_column',
            field=models.CharField(max_length=50, null=True),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='mapping_file',
            field=models.FileField(null=True, upload_to=''),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='quotechar',
            field=models.CharField(default='"', max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='platemapping',
            name='to_column',
            field=models.CharField(max_length=50, null=True),
        ),
    ]
