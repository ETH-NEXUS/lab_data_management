# Generated by Django 4.1.7 on 2023-03-30 08:44

from django.db import migrations


class Migration(migrations.Migration):
    initial = True

    def readFromFile(file: str) -> str:
        with open(file, "r") as f:
            return f.read()

    operations = [migrations.RunSQL(readFromFile("./pg_ext/db_scripts/functions.sql"))]
