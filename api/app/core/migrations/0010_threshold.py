# Generated by Django 4.2.3 on 2024-03-13 11:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0009_wellwithdrawal_current_amount_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='Threshold',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('dmso', models.FloatField(default=80)),
                ('amount', models.FloatField(default=2.5)),
            ],
        ),
    ]