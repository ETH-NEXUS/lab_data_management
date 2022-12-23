# Generated by Django 4.1.4 on 2022-12-20 08:56

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0023_alter_well_compounds'),
    ]

    operations = [
        migrations.AlterField(
            model_name='project',
            name='name',
            field=models.CharField(max_length=50, unique=True),
        ),
        migrations.AlterUniqueTogether(
            name='experiment',
            unique_together={('name', 'project')},
        ),
    ]
