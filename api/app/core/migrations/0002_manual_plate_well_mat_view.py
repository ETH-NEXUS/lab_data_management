# Generated by Django 4.1.7 on 2023-04-03 09:37

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [("pg_ext", "0001_manual_add_functions"), ("core", "0001_initial")]

    def readFromFile(file: str) -> str:
        with open(file, "r") as f:
            return f.read()

    operations = [
        migrations.RunSQL(readFromFile("./core/db_scripts/plate_well_mat_views.sql"))
    ]