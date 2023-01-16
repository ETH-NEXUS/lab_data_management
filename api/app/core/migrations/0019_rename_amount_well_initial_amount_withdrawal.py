# Generated by Django 4.1.4 on 2022-12-19 13:46

import core.models
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0018_remove_plate_project_experiment_plate_experiment'),
    ]

    operations = [
        migrations.RenameField(
            model_name='well',
            old_name='amount',
            new_name='initial_amount',
        ),
        migrations.CreateModel(
            name='Withdrawal',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('amount', models.FloatField()),
                ('well', models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='withdrawals', to='core.well')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
