from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [("core", "0017_plate_use_as_template_to_select")]

    def readFromFile(file: str) -> str:
        with open(file, "r") as f:
            return f.read()

    operations = [migrations.RunSQL(readFromFile("./core/db_scripts/functions.sql"))]
