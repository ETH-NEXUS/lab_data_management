from django.db import migrations

# Import looks a well type up in capitals (WellType.by_name), so "Nref" in a file is "NREF"
NAME = "NREF"
DESCRIPTION = "Negative Control Reference"


def add_nref(apps, schema_editor):
    WellType = apps.get_model("core", "WellType")
    # A new database gets its well types, NREF too, from the fixture well_types
    # after the migrations (`db init`); adding it here would take the id of "C"
    if not WellType.objects.exists():
        return
    WellType.objects.get_or_create(name=NAME, defaults={"description": DESCRIPTION})


class Migration(migrations.Migration):

    dependencies = [
        ("core", "0022_delete_location"),
    ]

    operations = [
        # Going back keeps NREF: wells may use it already
        migrations.RunPython(add_nref, migrations.RunPython.noop),
    ]
