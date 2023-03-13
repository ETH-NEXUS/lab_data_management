# Generated by Django 4.1.7 on 2023-03-13 08:52

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0007_alter_barcodespecification_experiment_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='barcodespecification',
            name='experiment',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='barcode_specifications', to='core.experiment'),
        ),
    ]
