# Generated by Django 4.2.3 on 2024-03-12 16:04

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compoundlib', '0014_remove_compound_identifier_remove_compound_library'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compound',
            name='structure',
            field=models.TextField(blank=True, null=True),
        ),
    ]