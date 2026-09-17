"""
Tests for the small value conversions used by the mappers.
"""

from django.test import SimpleTestCase

from datetime import datetime

from importer.mappers.values import convert_sci_to_float, parse_c10_datetime


class ParseC10DatetimeTest(SimpleTestCase):
    def test_the_date_and_time_of_a_txt_file(self):
        self.assertEqual(
            datetime(2024, 10, 14, 12, 45, 28),
            parse_c10_datetime("10/14/2024", "12:45:28"),
        )
        self.assertEqual(
            datetime(2025, 1, 5, 9, 5, 3), parse_c10_datetime("1/5/2025", "9:05:03")
        )

    def test_the_date_and_time_of_a_file_name(self):
        expected = datetime(2024, 10, 14, 12, 54, 55)
        self.assertEqual(expected, parse_c10_datetime("241014", "125455"))
        self.assertEqual(expected, parse_c10_datetime("20241014", "125455"))

    def test_an_unknown_format_gives_none(self):
        self.assertIsNone(parse_c10_datetime("14.10.2024", "12:45:28"))
        self.assertIsNone(parse_c10_datetime("2024101", "125455"))
        self.assertIsNone(parse_c10_datetime("10/14/2024", "12:45:28 PM"))

    def test_an_impossible_date_gives_none(self):
        self.assertIsNone(parse_c10_datetime("14/10/2024", "12:45:28"))
        self.assertIsNone(parse_c10_datetime("241314", "125455"))


class ConvertSciToFloatTest(SimpleTestCase):
    def test_numbers_are_converted(self):
        self.assertEqual(1500.0, convert_sci_to_float("1.5E+03"))
        self.assertEqual(409.0, convert_sci_to_float("409"))

    def test_text_gives_none(self):
        self.assertIsNone(convert_sci_to_float("abc"))
