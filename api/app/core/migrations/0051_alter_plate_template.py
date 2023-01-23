# Generated by Django 4.1.4 on 2023-01-16 13:33

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('platetemplate', '0002_alter_platetemplate_name'),
        ('core', '0050_remove_plate_check_only_library_or_experiment_or_template_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='plate',
            name='template',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='plates', to='platetemplate.platetemplate'),
        ),
    ]
