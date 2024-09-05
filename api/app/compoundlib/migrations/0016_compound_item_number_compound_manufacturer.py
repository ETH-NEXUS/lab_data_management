# Generated by Django 4.2.3 on 2024-09-05 11:44

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compoundlib', '0015_alter_compound_structure'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='item_number',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='manufacturer',
            field=models.TextField(blank=True, null=True),
        ),
    ]
