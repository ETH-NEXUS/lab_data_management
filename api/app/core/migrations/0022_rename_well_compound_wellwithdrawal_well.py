# Generated by Django 4.1.4 on 2022-12-19 15:08

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0021_wellwithdrawal_and_more'),
    ]

    operations = [
        migrations.RenameField(
            model_name='wellwithdrawal',
            old_name='well_compound',
            new_name='well',
        ),
    ]
