# Generated by Django 4.2.3 on 2024-03-14 09:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('problems', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='problem',
            name='show',
            field=models.BooleanField(default=True),
        ),
    ]