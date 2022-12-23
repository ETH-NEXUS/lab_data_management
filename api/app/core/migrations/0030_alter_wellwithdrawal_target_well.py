# Generated by Django 4.1.4 on 2022-12-22 13:01

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0029_alter_wellwithdrawal_target_well'),
    ]

    operations = [
        migrations.AlterField(
            model_name='wellwithdrawal',
            name='target_well',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='donors', to='core.well'),
        ),
    ]
