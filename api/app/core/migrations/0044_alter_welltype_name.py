# Generated by Django 4.1.4 on 2023-01-16 09:57

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0043_welltype_well_type'),
    ]

    operations = [
        migrations.AlterField(
            model_name='welltype',
            name='name',
            field=models.CharField(db_index=True, max_length=50),
        ),
    ]
