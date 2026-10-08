from django.db import migrations

# The views are created again from the same file as in 0002 and 0024. The fix in
# the file: the time points of a plate come from its overall stats, so they are
# in order also when well types have different reads. Running after pg_ext 0003
# also fills the views with the json_merge that does not search whole lists.
VIEWS_FILE = "./core/db_scripts/plate_well_mat_views.sql"


def read_views_file() -> str:
    with open(VIEWS_FILE, "r") as file:
        return file.read()


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0024_experiment_reads_counted_per_label"),
        ("pg_ext", "0003_json_merge_without_list_search"),
    ]

    operations = [
        migrations.RunSQL(read_views_file(), reverse_sql=migrations.RunSQL.noop),
    ]
