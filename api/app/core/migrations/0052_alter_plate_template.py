# Generated by Django 4.1.4 on 2023-01-16 14:19

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('platetemplate', '0003_alter_platetemplatecategory_name_and_more'),
        ('core', '0051_alter_plate_template'),
    ]

    operations = [
        migrations.AlterField(
            model_name='plate',
            name='template',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='plate', to='platetemplate.platetemplate'),
        ),
    ]
