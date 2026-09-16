"""
Tests for reading C10 reader/imager files (.txt and .xlsx) into measurement data.
"""

import shutil
import tempfile
from os.path import join

from django.test import SimpleTestCase
from openpyxl import Workbook

from importer.mappers import MicroscopeMapper

# A shortened real C10 txt export (the lines end with CRLF).
TXT_LINES = [
    "",
    "",
    "Software Version\t3.14.03",
    "",
    "Experiment File Path:\tG:\\Reader_ExperimentFiles\\Pruschy_CellTiterGlo3D_241014_125455.xpt",
    "Plate Number\tPlate 1",
    "Date\t10/14/2024",
    "Time\t12:45:28",
    "Reader Type:\tCytationC10",
    "Procedure Details",
    "Read\tLuminescence Endpoint",
    "\tFull Plate",
    "\tIntegration Time: 0:01.00 (MM:SS.ss)",
    "Automatic gain values\t ",
    "Gain(Lum)\t214",
    "Actual Temperature:\t28.6",
    "",
    "Results",
    "",
    "Well\tLum",
    "A1\t16727",
    "A2\t1.5E+03",
    "P24\t409",
    "",
]

# A made-up xlsx with the parts the parser looks for: metadata at the top,
# a "Layout" block with POS/NEG wells and a "Results" table.
XLSX_ROWS = [
    ["Plate Number", "Plate 1"],
    ["Date", "10/14/2024"],
    ["Read", "Luminescence Endpoint", "Gain: 214"],
    ["Layout"],
    [None, 1, 2],
    [None, "A", "POS", "SMP"],
    [None, "B", "NEG", "SMP"],
    [],
    ["Results"],
    [None, "Well ID", "Well", "Lum"],
    [None, "SPL1", "A1", 16727],
    [None, "SPL2", "B1", "1.5E+03"],
]


class MicroscopeParseTest(SimpleTestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def write_txt(self):
        path = join(self.folder, "241014_125455_241008MP-1_1.txt")
        with open(path, "w", newline="") as file:
            file.write("\r\n".join(TXT_LINES))
        return path

    def write_xlsx(self):
        path = join(self.folder, "241014_125455_241008MP-1_2.xlsx")
        workbook = Workbook()
        for row in XLSX_ROWS:
            workbook.active.append(row)
        workbook.save(path)
        return path

    def test_a_txt_file_is_read(self):
        data = MicroscopeMapper().parse(self.write_txt())

        self.assertEqual("241008MP-1_1", data["barcode"])
        # Date and time come from the file content, not from the file name
        self.assertEqual("10/14/2024", data["date"])
        self.assertEqual("12:45:28", data["time"])
        self.assertEqual({}, data["layout"])
        self.assertEqual(
            [
                {"Well": "A1", "Lum": "16727"},
                {"Well": "A2", "Lum": "1.5E+03"},
                {"Well": "P24", "Lum": "409"},
            ],
            data["results"],
        )
        self.assertEqual(
            {
                "Software Version": "3.14.03",
                "Experiment File Path:": "G:\\Reader_ExperimentFiles\\Pruschy_CellTiterGlo3D_241014_125455.xpt",
                "Plate Number": "Plate 1",
                "Date": "10/14/2024",
                "Time": "12:45:28",
                "Reader Type:": "CytationC10",
                "Read": "Luminescence Endpoint",
                "Integration Time": "0:01.00 (MM:SS.ss)",
                "Gain(Lum)": "214",
                "Actual Temperature:": "28.6",
            },
            data["metadata"],
        )

    def test_a_txt_file_uses_the_given_measurement_name(self):
        data = MicroscopeMapper().parse(self.write_txt(), measurement_name="Label1")

        self.assertEqual({"Well": "A1", "Label1": "16727"}, data["results"][0])

    def test_a_txt_file_with_measurement_name_none_labels_the_values_none(self):
        # Current behavior: `map` passes measurement_name=None when no name is
        # given, and then the values are stored under the key None.
        data = MicroscopeMapper().parse(self.write_txt(), measurement_name=None)

        self.assertEqual({"Well": "A1", None: "16727"}, data["results"][0])

    def test_an_xlsx_file_is_read(self):
        data = MicroscopeMapper().parse(self.write_xlsx())

        self.assertEqual("241008MP-1_2", data["barcode"])
        # Date and time come from the file name
        self.assertEqual("241014", data["date"])
        self.assertEqual("125455", data["time"])
        self.assertEqual({"A1": "P", "B1": "N"}, data["layout"])
        self.assertEqual(
            [
                {"Well ID": "SPL1", "Well": "A1", "Lum": 16727},
                {"Well ID": "SPL2", "Well": "B1", "Lum": "1.5E+03"},
            ],
            data["results"],
        )
        self.assertEqual(
            {
                "Plate Number": ["Plate 1"],
                "Date": ["10/14/2024"],
                "Read": ["Luminescence Endpoint"],
                "Gain": "214",
                "Layout": ["1", "2", "A", "POS", "SMP", "B", "NEG", "SMP"],
                "Results": [
                    "Well ID",
                    "Well",
                    "Lum",
                    "SPL1",
                    "A1",
                    "16727",
                    "SPL2",
                    "B1",
                    "1.5E+03",
                ],
            },
            data["metadata"],
        )

    def test_a_file_name_without_date_and_time_is_refused(self):
        path = join(self.folder, "241008MP-1_1.txt")
        with open(path, "w") as file:
            file.write("Results\n")

        with self.assertRaises(ValueError):
            MicroscopeMapper().parse(path)
