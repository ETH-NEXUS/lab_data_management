# Generated by Django 4.1.7 on 2023-03-17 11:58

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0009_platemapping_evaluation'),
    ]

    operations = [
        migrations.AlterField(
            model_name='measurement',
            name='created_at',
            field=models.DateTimeField(auto_now_add=True),
        ),
        migrations.CreateModel(
            name='MeasurementAssignment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('assignment_timestamp', models.DateTimeField(auto_now_add=True)),
                ('status', models.CharField(default='pending', max_length=50)),
                ('filename', models.TextField()),
                ('measurement_file', models.FileField(null=True, upload_to='')),
                ('measurement_timestamp', models.DateTimeField(null=True)),
                ('measurement', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='assignments', to='core.measurement')),
                ('metadata', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='assignments', to='core.measurementmetadata')),
                ('plate', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='assignments', to='core.plate')),
            ],
        ),
    ]
