"""
Tests for the small value conversions used by the mappers.
"""

from django.test import SimpleTestCase

from importer.mappers import values
from importer.mappers.values import convert_sci_to_float, convert_string_to_datetime


class ConvertStringToDatetimeTest(SimpleTestCase):
    def test_a_date_and_time_in_yyyymmdd_hhmmss_format_is_converted(self):
        self.assertEqual(
            "2024-10-14T12:54:55+00:00",
            convert_string_to_datetime("20241014", "125455"),
        )

    def test_other_formats_fall_back_to_the_time_the_module_was_loaded(self):
        # Current behavior: the C10 dates "241014" (file name) and
        # "10/14/2024" (txt content) do not fit, so GLOBAL_NOW is used.
        fallback = values.GLOBAL_NOW.isoformat()

        self.assertEqual(fallback, convert_string_to_datetime("241014", "125455"))
        self.assertEqual(fallback, convert_string_to_datetime("10/14/2024", "12:45:28"))


class ConvertSciToFloatTest(SimpleTestCase):
    def test_numbers_are_converted(self):
        self.assertEqual(1500.0, convert_sci_to_float("1.5E+03"))
        self.assertEqual(409.0, convert_sci_to_float("409"))

    def test_text_gives_none(self):
        self.assertIsNone(convert_sci_to_float("abc"))
