"""
The plate barcode and the time of the measurement from the name of an M1000
file, in the format chosen on the management page.
"""

from datetime import datetime

from django.core.management.base import CommandError
from django.test import SimpleTestCase

from importer.mappers.base import SkipFile
from importer.mappers.m1000_file_names import read_file_name


class ReadFileNameTest(SimpleTestCase):
    def test_the_older_format_with_date_and_time(self):
        self.assertEqual(
            {"barcode": "demo_1", "measured_at": datetime(2024, 6, 10, 12, 12, 12)},
            read_file_name("/data/20240610-121212_demo_1.asc", "date_time_barcode"),
        )

    def test_the_older_format_without_date(self):
        self.assertEqual(
            {"barcode": "demo_1", "measured_at": None},
            read_file_name("/data/demo_1.asc", "date_time_barcode"),
        )

    def test_the_older_format_with_the_6_digit_date_of_the_evaluated_files(self):
        # As the lab names them: 09/30/26 15:46:54, month first. Read with the
        # 8 digit pattern it was the year 930
        self.assertEqual(
            {
                "barcode": "RKS_300926_1",
                "measured_at": datetime(2026, 9, 30, 15, 46, 54),
            },
            read_file_name("/data/093026-154654_RKS_300926_1.asc", "date_time_barcode"),
        )

    def test_the_older_format_with_another_number_of_digits_is_refused(self):
        with self.assertRaisesMessage(
            CommandError,
            "The date 93026 in the file name 93026-154654_RKS_1.asc has 5 digits, "
            "but 8 or 6 are expected",
        ):
            read_file_name("93026-154654_RKS_1.asc", "date_time_barcode")

    def test_without_a_format_the_older_one_is_used(self):
        self.assertEqual(
            "240716MP-1_2",
            read_file_name("20240722-125833_240716MP-1_2.asc", None)["barcode"],
        )

    def test_the_newer_format_as_the_lab_names_its_files(self):
        # The date is month, day, year: 09/30/26
        self.assertEqual(
            {
                "barcode": "RKS_300926_3",
                "measured_at": datetime(2026, 9, 30, 16, 54, 54),
            },
            read_file_name("/data/RKS_300926_3_093026_165454.asc", "barcode_date_time"),
        )

    def test_another_name_in_the_newer_format_is_skipped(self):
        # The lab's folder has "30092026-001.asc" next to the files of the plates
        with self.assertRaisesMessage(
            SkipFile,
            "its name does not match the chosen format barcode_date_time, "
            "e.g. RKS_300926_3_093026_165454.asc.",
        ):
            read_file_name("30092026-001.asc", "barcode_date_time")

    def test_a_newer_name_in_the_older_format_is_refused(self):
        # Without this check the whole name would become the barcode of a new plate
        with self.assertRaisesMessage(
            CommandError,
            "The file name RKS_300926_2_093026_162008.asc looks like the format "
            "barcode_date_time, e.g. RKS_300926_3_093026_165454.asc. "
            "Choose that format.",
        ):
            read_file_name("RKS_300926_2_093026_162008.asc", "date_time_barcode")

    def test_an_unknown_format_is_refused(self):
        with self.assertRaisesMessage(CommandError, "Unknown file name format 'x'"):
            read_file_name("demo_1.asc", "x")

    def test_a_date_that_does_not_exist_is_no_date(self):
        # Month 13
        self.assertIsNone(
            read_file_name("RKS_300926_3_133026_165454.asc", "barcode_date_time")[
                "measured_at"
            ]
        )
