"""
Tests for the SDF amounts in nanoliter.
"""

from django.test import SimpleTestCase

from importer.sdf_amounts import amount_in_nanoliter, is_volume_column


class SdfAmountsTest(SimpleTestCase):
    def test_a_volume_in_microliter_is_stored_in_nanoliter(self):
        self.assertEqual(6000.0, amount_in_nanoliter("6"))
        self.assertEqual(24000.0, amount_in_nanoliter("24.0"))
        self.assertEqual(6000.0, amount_in_nanoliter(" 6 "))

    def test_an_empty_value_is_zero(self):
        self.assertEqual(0.0, amount_in_nanoliter(""))
        self.assertEqual(0.0, amount_in_nanoliter(None))
        self.assertEqual(0.0, amount_in_nanoliter(float("nan")))

    def test_a_value_that_is_not_a_number_gives_none(self):
        self.assertIsNone(amount_in_nanoliter("<24"))

    def test_only_the_vol_copy_columns_are_volumes(self):
        self.assertTrue(is_volume_column("Vol_Copy1"))
        self.assertFalse(is_volume_column("PLATE_AMOUNT1"))
