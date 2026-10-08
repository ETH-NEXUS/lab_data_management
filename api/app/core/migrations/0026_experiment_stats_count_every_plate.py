from django.db import migrations

# The views are created again from the same file as in 0002, 0024 and 0025. The fix
# in the file: the overall stats of an experiment take the values of every plate,
# also when two plates have exactly the same values in a read.
VIEWS_FILE = "./core/db_scripts/plate_well_mat_views.sql"


def read_views_file() -> str:
    with open(VIEWS_FILE, "r") as file:
        return file.read()


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0025_plate_timestamps_from_overall_stats"),
    ]

    operations = [
        migrations.RunSQL(read_views_file(), reverse_sql=migrations.RunSQL.noop),
    ]
