"""
json_merge (pg_ext) merges the values of the same key in a JSON object, e.g. the
timestamps of a label that the views list once per well type.
"""

from django.db import connection
from django.test import TestCase


def json_merge(text: str) -> dict:
    with connection.cursor() as cursor:
        cursor.execute("SELECT json_merge(%s::json)", [text])
        # The database driver reads the JSON result as a dict already
        return cursor.fetchone()[0]


class JsonMergeTest(TestCase):
    def test_lists_keep_their_order_and_lose_duplicates(self):
        merged = json_merge('{"Lum": ["10:00", "12:00"], "Lum": ["12:00", "14:00"]}')

        self.assertEqual({"Lum": ["10:00", "12:00", "14:00"]}, merged)

    def test_lists_of_objects_are_merged_too(self):
        merged = json_merge('{"a": [{"x": 1}], "a": [{"x": 1}, {"x": 2}]}')

        self.assertEqual({"a": [{"x": 1}, {"x": 2}]}, merged)

    def test_objects_are_merged_by_their_keys(self):
        merged = json_merge('{"Lum": {"C": 1}, "Lum": {"N1": 2}}')

        self.assertEqual({"Lum": {"C": 1, "N1": 2}}, merged)
