from django.db import migrations

# The views are dropped and created again from the same file as in 0002, so a new
# database and an old one end up with the same views. The fix in the file: the
# overall stats of an experiment count the reads of a plate per label. Creating
# them again also fills them with the fixed json_merge (pg_ext 0002).
VIEWS_FILE = "./core/db_scripts/plate_well_mat_views.sql"


def read_views_file() -> str:
    with open(VIEWS_FILE, "r") as file:
        return file.read()


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0023_add_well_type_nref"),
        # The views are filled with the json_merge that keeps the order of the timestamps
        ("pg_ext", "0002_json_merge_keeps_order"),
    ]

    operations = [
        migrations.RunSQL(read_views_file(), reverse_sql=migrations.RunSQL.noop),
    ]
