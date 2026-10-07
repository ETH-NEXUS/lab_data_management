"""
The background of the reference wells is subtracted from the other wells,
the reference wells are left out.
"""

from django.test import SimpleTestCase

from background_correction.calculation import subtract_background


class SubtractBackgroundTest(SimpleTestCase):
    def test_median_of_the_reference_wells_is_subtracted(self):
        values = {11: 10.0, 12: 4.0, 13: 2.0, 14: 6.0, 15: 100.0}

        corrected = subtract_background(values, {13, 14, 15}, "median")

        # Median of 2, 6 and 100 is 6
        self.assertEqual({11: 4.0, 12: -2.0}, corrected)

    def test_mean_of_the_reference_wells_is_subtracted(self):
        values = {11: 10.0, 12: 4.0, 13: 2.0, 14: 6.0, 15: 100.0}

        corrected = subtract_background(values, {13, 14, 15}, "mean")

        # Mean of 2, 6 and 100 is 36
        self.assertEqual({11: -26.0, 12: -32.0}, corrected)

    def test_reference_wells_without_a_value_are_ignored(self):
        values = {11: 10.0, 13: 2.0}

        corrected = subtract_background(values, {13, 14}, "median")

        self.assertEqual({11: 8.0}, corrected)
