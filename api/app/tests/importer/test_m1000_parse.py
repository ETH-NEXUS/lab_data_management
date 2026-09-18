"""
Tests for reading M1000 reader files (.asc) into measurement entries.
"""

import shutil
import tempfile
from datetime import datetime
from io import StringIO
from unittest import mock
from os.path import join

from django.core.management.base import CommandError
from django.test import SimpleTestCase

from importer.mappers import M1000Mapper

# A shortened demo file. The line endings are as in the real file: the value
# lines end with LF, the last value line and the footer with CRLF.
ASC_FILE_CONTENT = (
    b"A1\tSM1_1\t15\n"
    b"A2\tSM1_2\t999\n"
    b"B1\tSM1_25\t1.5E+02\n"
    b"P24\tSM1_384\t9\r\n"
    b"Date of measurement: 2011-11-11/Time of measurement: 11:11:11\r\n"
    b"test.mth\r\n"
    b"\r\n"
    b"infinite M1000 PRO\r\n"
    b"Instrument serial number: 1205002707\r\n"
    b"Plate\r\n"
    b"Plate Description: test\r\n"
    b"Plate with Cover: No\r\n"
    b"Barcode: No\r\n"
    b"  Part of Plate\r\n"
    b"  Range: A1:P24\r\n"
    b"    Demo\r\n"
    b"    Attenuation: NONE\r\n"
    b"    Integration time: 1000 ms\r\n"
    b"    Label: Label1\r\n"
    b"    Settle time: 0 ms\r\n"
    b"Meas. temperature: Raw data: 27.3 M-0C\r\n"
    b"Date: 2024-06-10, Time: 11:11:11\r\n"
)


class M1000ParseTest(SimpleTestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.path = join(self.folder, "20240610-121212_demo_1.asc")
        with open(self.path, "wb") as file:
            file.write(ASC_FILE_CONTENT)

    def tearDown(self):
        shutil.rmtree(self.folder)

    def parse(self):
        # Opened in text mode, the same way BaseMapper.run opens it
        with open(self.path, "r", encoding="ascii") as file:
            return M1000Mapper().parse(file, room_name=None)

    def test_the_value_lines_become_entries(self):
        results = self.parse()["entries"]

        self.assertEqual(
            [
                {"position": "A1", "identifier": "SM1_1", "values": [15.0]},
                {"position": "A2", "identifier": "SM1_2", "values": [999.0]},
                {"position": "B1", "identifier": "SM1_25", "values": [150.0]},
                {"position": "P24", "identifier": "SM1_384", "values": [9.0]},
            ],
            results,
        )

    def test_the_footer_and_the_file_name_give_the_plate_information(self):
        data = self.parse()
        del data["entries"]

        self.assertEqual(
            {
                "barcode": "demo_1",
                "measurement_date": datetime(2011, 11, 11, 11, 11, 11),
                "plate_description": "test",
                "meta_data": [
                    {
                        "Attenuation": "NONE",
                        "Integration time": "1000 ms",
                        "Label": "Label1",
                        "Settle time": "0 ms",
                    }
                ],
            },
            data,
        )

    def test_a_value_that_only_starts_like_a_number_stops_the_file(self):
        # "12abc" matches the number pattern at its start, but is no number
        with self.assertRaisesMessage(
            CommandError, "The value '12abc' of well A1 is not a number."
        ):
            M1000Mapper().read_value_line(["A1", "SM1_1", "12abc"], 0, 1)

    def test_a_column_that_is_no_number_keeps_its_place(self):
        # "OVER" means the reader saw more light than it can measure
        entry = M1000Mapper().read_value_line(["A1", "SM1_1", "OVER", "15"], 0, 1)

        self.assertEqual([None, 15.0], entry["values"])


class M1000DetermineIndexesTest(SimpleTestCase):
    def determine(self, text):
        """(position column, identifier column) and where the file is afterwards."""
        file = StringIO(text)
        file.name = "test.asc"
        indexes = M1000Mapper().determine_indexes(file)
        return indexes, file.tell()

    def test_the_position_column_comes_first(self):
        self.assertEqual(((0, 1), 0), self.determine("A1\tSM1_1\t15\nA2\tSM1_2\t16\n"))

    def test_the_identifier_column_comes_first(self):
        # "NC1" looks like a position too, so the first line is ambiguous and
        # the second line decides.
        text = "NC1\tA1\t11115\nSM1_1\tA3\t12154\n"

        self.assertEqual(((1, 0), 0), self.determine(text))

    def test_a_header_line_above_the_values_is_skipped(self):
        text = "Well positions\tLayout\tAcceptor\tDonor\t\nA1\tSM1_1\t24672\t7395\t\n"

        self.assertEqual(((0, 1), 0), self.determine(text))

    def test_spaces_around_the_values_are_ignored(self):
        self.assertEqual(((0, 1), 0), self.determine("  A1   SM1_1  15  \n"))

    def test_a_file_without_an_unambiguous_line_is_refused(self):
        file = StringIO("A1\t15\n\nB1\tSM1_2\tSM1_3\t4\n")
        file.name = "test.asc"

        with self.assertRaises(CommandError) as raised:
            M1000Mapper().determine_indexes(file)

        self.assertEqual(
            "File has not the desired format: test.asc", str(raised.exception)
        )
        self.assertEqual(0, file.tell())


class M1000MeasurementDateTest(SimpleTestCase):
    def parse(self, name, footer=""):
        file = StringIO("A1\tSM1_1\t15\n" + footer)
        file.name = name
        return M1000Mapper().parse(file)

    def test_the_date_of_the_footer_is_used(self):
        data = self.parse(
            "/data/20240610-121212_demo_1.asc",
            "Date of measurement: 2024-06-11/Time of measurement: 08:30:00\n",
        )

        self.assertEqual(datetime(2024, 6, 11, 8, 30, 0), data["measurement_date"])

    def test_without_a_date_in_the_footer_the_file_name_is_used(self):
        data = self.parse("/data/20240610-121212_demo_1.asc")

        self.assertEqual(datetime(2024, 6, 10, 12, 12, 12), data["measurement_date"])

    @mock.patch("importer.mappers.m1000.message")
    def test_without_any_date_the_time_of_the_import_is_used(self, message):
        data = self.parse("/data/demo_1.asc")

        self.assertIsNotNone(data["measurement_date"])
        self.assertIn(
            "Neither the file nor its name says when it was measured",
            message.call_args.args[0],
        )
        self.assertEqual("warning", message.call_args.args[1])
