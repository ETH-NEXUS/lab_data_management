"""
Tests for reading and writing well positions in the "A3" notation.
"""

from django.test import SimpleTestCase

from core.utils.plates.positions import PositionMapper, PositionMappingError


class PositionMapperTest(SimpleTestCase):
    def test_a_position_is_turned_into_row_and_column(self):
        self.assertEqual((1, 3), PositionMapper.map("A3"))
        self.assertEqual((2, 12), PositionMapper.map("b12"))
        self.assertEqual((27, 1), PositionMapper.map("AA1"))

    def test_row_and_column_are_turned_into_a_position(self):
        self.assertEqual("A3", PositionMapper.unmap(1, 3))
        self.assertEqual("AA1", PositionMapper.unmap(27, 1))

    def test_an_invalid_position_raises_a_position_mapping_error(self):
        # This used to end in a TypeError, because the error class did not
        # inherit from Exception.
        with self.assertRaises(PositionMappingError) as raised:
            PositionMapper.map("??")
        self.assertEqual(
            "Cannot convert position to row, col: ??", raised.exception.message
        )
        self.assertEqual("Cannot convert position to row, col: ??", str(raised.exception))
