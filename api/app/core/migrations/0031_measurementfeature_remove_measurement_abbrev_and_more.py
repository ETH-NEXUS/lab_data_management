# Generated by Django 4.1.4 on 2022-12-23 09:02

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0030_alter_wellwithdrawal_target_well'),
    ]

    operations = [
        migrations.CreateModel(
            name='MeasurementFeature',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=50, null=True, verbose_name='measurement')),
                ('abbrev', models.CharField(blank=True, max_length=4, null=True)),
                ('unit', models.CharField(blank=True, max_length=10, null=True)),
            ],
        ),
        migrations.RemoveField(
            model_name='measurement',
            name='abbrev',
        ),
        migrations.RemoveField(
            model_name='measurement',
            name='name',
        ),
        migrations.RemoveField(
            model_name='measurement',
            name='unit',
        ),
        migrations.AddField(
            model_name='measurement',
            name='feature',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.RESTRICT, related_name='measurements', to='core.measurementfeature'),
            preserve_default=False,
        ),
    ]
