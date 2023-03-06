from copy import deepcopy
from os import makedirs
from os.path import join
from pathlib import Path

from django.test import TestCase

from core.models import PlateDimension
from .helper import sameSchema, row_col_from_wells, closest, row_col_from_name
from .mappers import BaseMapper

from pathlib import Path
from core.helper import posToAlphaChar
from core.models import PlateDimension
from .mappers import BaseMapper, EchoMapper


class HelperTest(TestCase):
    def test_sameSchemaTrue(self):
        a = {"a": {"b": {"c": "xxx"}, "x": {"y": 111}}}
        a = {"a": {"b": {"c": "xxx"}, "x": {"y": 111}}}
        b = deepcopy(a)
        self.assertTrue(sameSchema(a, b))

    def test_sameSchemaFalse(self):
        a = {"a": {"b": {"c": "xxx"}, "x": {"y": 111}}}
        a = {"a": {"b": {"c": "xxx"}, "x": {"y": 111}}}
        b = deepcopy(a)
        del b["a"]["x"]["y"]
        b["a"]["x"]["z"] = "changed"
        del b["a"]["x"]["y"]
        b["a"]["x"]["z"] = "changed"
        self.assertFalse(sameSchema(a, b))

    def test_row_col_from_wells(self):
        """Test row and col from wells"""
        """Test row and col from wells"""
        r, c = row_col_from_wells(96)
        self.assertEqual(8, r)
        self.assertEqual(12, c)
        r, c = row_col_from_wells(384)
        self.assertEqual(16, r)
        self.assertEqual(24, c)
        r, c = row_col_from_wells(1536)
        self.assertEqual(32, r)
        self.assertEqual(48, c)

    def test_closest(self):
        """Test the closest function"""
        """Test the closest function"""
        lst = (12, 24, 48)
        self.assertEqual(12, closest(6, lst))
        self.assertEqual(12, closest(18, lst))
        self.assertEqual(24, closest(19, lst))
        self.assertEqual(24, closest(36, lst))
        self.assertEqual(48, closest(37, lst))
        self.assertEqual(24, closest(24, lst))
        self.assertEqual(48, closest(100, lst))

    def test_row_col_from_name(self):
        r, c = row_col_from_name("Greiner_384PS_781090")
        self.assertEqual(16, r)
        self.assertEqual(24, c)

    def test_row_col_from_name_exception(self):
        try:
            row_col_from_name("Greiner_8989PS_781090")
            self.fail("Should raise ValueError")
        except ValueError as ve:
            self.assertEqual(
                "Cannot determine plate dimension from name: Greiner_8989PS_781090.",
                str(ve),
            )


class ConvertPositionToIndexTests(TestCase):
    fixtures = ["plate_dimensions"]
    fixtures = ["plate_dimensions"]

    def setUp(self):
        self.dimension_96 = PlateDimension.objects.get(name="dim_96_8x12")
        self.dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        self.dimension_1536 = PlateDimension.objects.get(name="dim_1536_32x48")
        self.dimension_96 = PlateDimension.objects.get(name="dim_96_8x12")
        self.dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        self.dimension_1536 = PlateDimension.objects.get(name="dim_1536_32x48")

    def test_single_letter(self):
        """
        """
        Test a position with a single-letter row label
        """
        """

        position = "B3"
        expected_index = 26
        self.assertEqual(self.dimension_384.position(position), expected_index)
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_first_column(self):
        """
        """
        Test a position in the first column
        """
        """

        position = "A1"
        expected_index = 0
        self.assertEqual(self.dimension_96.position(position), expected_index)

    def test_random_column(self):
        """
        """
        Test a position in the last column
        """
        """

        position = "D3"
        expected_index = 74
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_last_column(self):
        """
        """
        Test a position in the last column
        """
        """

        position = "AA11"
        expected_index = 1258
        self.assertEqual(self.dimension_1536.position(position), expected_index)

    def test_index_to_letter(self):
        """
        """
        Test a position in the last column
        """
        """

        index = 27
        expected_letter = "AA"
        self.assertEqual(posToAlphaChar(index), expected_letter)
        expected_letter = "AA"
        self.assertEqual(posToAlphaChar(index), expected_letter)
        self.assertEqual(posToAlphaChar(28), "AB")
        self.assertEqual(posToAlphaChar(52), "AZ")
        self.assertEqual(posToAlphaChar(703), "AAA")
        self.assertEqual(posToAlphaChar(73116), "DDDD")


class MapperTests(TestCase):
    fixtures = ["plate_dimensions", "well_types", "test/compound_library"]
    fixtures = ["plate_dimensions", "well_types", "test/compound_library"]

    TEST_DATA_FOLDER = "./temp"
    ECHO_DIR = join(TEST_DATA_FOLDER, "echo")
    BARCODES = ["P1", "P2", "P3"]
    TEST_DATA_FOLDER = "./temp"
    ECHO_DIR = join(TEST_DATA_FOLDER, "echo")
    BARCODES = ["P1", "P2", "P3"]

    def create_echo_test_data(self):
        for barcode in self.BARCODES:
            _dir = join(self.ECHO_DIR, barcode)
            makedirs(_dir, exist_ok=True)
            Path(join(_dir, "something.csv")).touch()
            Path(join(_dir, "dfg-transfer-file.tsv")).touch()
            with open(join(_dir, "ID-123-transfer-Echo_01_123.csv"), "w") as ef:
                ef.write(
                    f"""
            Path(join(_dir, "something.csv")).touch()
            Path(join(_dir, "dfg-transfer-file.tsv")).touch()
            with open(join(_dir, "ID-123-transfer-Echo_01_123.csv"), "w") as ef:
                ef.write(
                    f"""
                    Run ID,14618
                    Run Date/Time,02/09/2021 10:30:32
                    Application Name,
                    Application Version,1.0.0.0
                    Protocol Name,
                    User Name,cellario


                    [DETAILS]
                    Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Source Concentration,Source Concentration Units,Destination Plate Name,Destination Plate Barcode,Destination Well,Destination Concentration,Destination Concentration Units,Compound Name,Transfer Volume,Actual Volume,Transfer Status,Current Fluid Height,Current Fluid Volume,% DMSO
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A10,0,N/A,Corning_384_3577,{barcode},A10,0,N/A,N/A,30,30,,2.003,9.959,99.203
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A11,0,N/A,Corning_384_3577,{barcode},A11,0,N/A,N/A,30,30,,1.968,9.786,99.008
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A12,0,N/A,Corning_384_3577,{barcode},A12,0,N/A,N/A,30,30,,2.001,9.959,99.348
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A13,0,N/A,Corning_384_3577,{barcode},A13,0,N/A,N/A,30,30,,2.013,9.983,98.519
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A14,0,N/A,Corning_384_3577,{barcode},A14,0,N/A,N/A,30,30,,2.001,9.973,99.701
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A15,0,N/A,Corning_384_3577,{barcode},A15,0,N/A,N/A,30,30,,1.964,9.793,99.643
                    """
                )
                    """
                )

    def setUp(self):
        self.create_echo_test_data()

    def tearDown(self):
        # TODO: Delete files
        pass

    def test_get_files(self):
        mapper = BaseMapper()
        self.assertListEqual(
            [
                join(self.ECHO_DIR, barcode, "ID-123-transfer-Echo_01_123.csv")
                for barcode in self.BARCODES
            ],
            mapper.get_files(join(self.ECHO_DIR, "**", "*-transfer-*.csv")),
        )
