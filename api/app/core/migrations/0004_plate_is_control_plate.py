# Generated by Django 4.2.2 on 2023-07-11 17:19

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0003_experimentdetail'),
    ]

    operations = [
        migrations.AddField(
            model_name='plate',
            name='is_control_plate',
            field=models.BooleanField(blank=True, default=False, null=True),
        ),
    ]
