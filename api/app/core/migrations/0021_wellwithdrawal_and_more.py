# Generated by Django 4.1.4 on 2022-12-19 14:29

import core.models
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0020_rename_from_well_well_source_wells_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='WellWithdrawal',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.models.AutoDateTimeField(editable=False)),
                ('amount', models.FloatField()),
                ('well_compound', models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='well_withdrawals', to='core.well')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RenameField(
            model_name='wellcompound',
            old_name='initial_amount',
            new_name='amount',
        ),
        migrations.DeleteModel(
            name='Withdrawal',
        ),
    ]
