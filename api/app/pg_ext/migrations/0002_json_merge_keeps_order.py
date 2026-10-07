from django.db import migrations

# Only the body of json_merge is replaced: the materialized views use it, and
# dropping it (as functions.sql does) would drop them. The same function as in
# functions.sql, where a new database gets it from.
JSON_MERGE_KEEPS_ORDER = """
CREATE OR REPLACE FUNCTION json_merge(jdata json) RETURNS json AS $$
import json

def merge_on_duplicate_keys(ordered_pairs):
        d = {}
        for k, v in ordered_pairs:
                if k in d:
                        if isinstance(d[k], dict):
                                d[k].update(v)
                        elif isinstance(d[k], list):
                                # Keeps the order (e.g. the timestamps of a label), a set would not
                                for item in v:
                                        if item not in d[k]:
                                                d[k].append(item)
                        else:
                                d[k] += f", {v}"
                else:
                        d[k] = v
        return d

return json.dumps(json.loads(jdata, object_pairs_hook=merge_on_duplicate_keys))
$$ LANGUAGE plpython3u;
"""


class Migration(migrations.Migration):

    dependencies = [
        ("pg_ext", "0001_manual_add_functions"),
    ]

    operations = [
        migrations.RunSQL(JSON_MERGE_KEEPS_ORDER, reverse_sql=migrations.RunSQL.noop),
    ]
