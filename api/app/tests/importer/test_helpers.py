from copy import deepcopy
from django.test import TestCase
from core.utils.plates.positions import posToAlphaChar
from core.models import PlateDimension
from importer.helper import sameSchema, row_col_from_wells, closest, row_col_from_name


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

    def setUp(self):
        self.dimension_96 = PlateDimension.objects.get(name="dim_96_8x12")
        self.dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        self.dimension_1536 = PlateDimension.objects.get(name="dim_1536_32x48")
        self.dimension_96 = PlateDimension.objects.get(name="dim_96_8x12")
        self.dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        self.dimension_1536 = PlateDimension.objects.get(name="dim_1536_32x48")

    def test_single_letter(self):
        """
        Test a position with a single-letter row label
        """

        position = "B3"
        expected_index = 26
        self.assertEqual(self.dimension_384.position(position), expected_index)
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_first_column(self):
        """
        Test a position in the first column
        """

        position = "A1"
        expected_index = 0
        self.assertEqual(self.dimension_96.position(position), expected_index)

    def test_random_column(self):
        """
        Test a position in the last column
        """

        position = "D3"
        expected_index = 74
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_last_column(self):
        """
        Test a position in the last column
        """

        position = "AA11"
        expected_index = 1258
        self.assertEqual(self.dimension_1536.position(position), expected_index)

    def test_index_to_letter(self):
        """
        Test a position in the last column
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
