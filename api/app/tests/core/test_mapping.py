from django.test import TestCase
from core.models import PlateDimension
from core.helper import charToAlphaPos


class MappingTest(TestCase):
    def test_charToAlphaPos(self):
        self.assertEqual(1, charToAlphaPos("A"))
        self.assertEqual(16, charToAlphaPos("P"))
        self.assertEqual(17, charToAlphaPos("q"))
        self.assertEqual(26, charToAlphaPos("z"))
        self.assertRaises(ValueError, lambda: charToAlphaPos("123"))

    def test_positionMapping(self):
        plateDimension = PlateDimension.objects.create(name="bla", rows=2, cols=3)
        pos = plateDimension.position("B2")
        self.assertEqual(4, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual("B2", hr_pos)

        plateDimension = PlateDimension.objects.create(name="bla", rows=4, cols=5)
        pos = plateDimension.position("C3")
        self.assertEqual(12, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual("C3", hr_pos)

        pos = plateDimension.position("D2")
        self.assertEqual(16, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual("D2", hr_pos)

        plateDimension = PlateDimension.objects.create(name="bla", rows=16, cols=22)
        pos = plateDimension.position("A03")
        self.assertEqual(2, pos)
        hr_pos = plateDimension.hr_position(pos)
        self.assertEqual("A3", hr_pos)
