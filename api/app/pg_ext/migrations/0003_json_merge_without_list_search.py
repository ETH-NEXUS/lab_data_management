from django.db import migrations

# json_merge without searching the whole list for every item: with many time points
# the refresh of the views slowed down. Only the body is replaced, as in 0002.
JSON_MERGE_WITHOUT_LIST_SEARCH = """
CREATE OR REPLACE FUNCTION json_merge(jdata json) RETURNS json AS $$
import json

def merge_on_duplicate_keys(ordered_pairs):
        d = {}
        for k, v in ordered_pairs:
                if k in d:
                        if isinstance(d[k], dict):
                                d[k].update(v)
                        elif isinstance(d[k], list):
                                # Keeps the order (e.g. the timestamps of a label), a set alone would not.
                                # The set of the items as text finds a duplicate without searching the list.
                                seen = {json.dumps(item, sort_keys=True) for item in d[k]}
                                for item in v:
                                        key = json.dumps(item, sort_keys=True)
                                        if key not in seen:
                                                seen.add(key)
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
        ("pg_ext", "0002_json_merge_keeps_order"),
    ]

    operations = [
        migrations.RunSQL(
            JSON_MERGE_WITHOUT_LIST_SEARCH, reverse_sql=migrations.RunSQL.noop
        ),
    ]
