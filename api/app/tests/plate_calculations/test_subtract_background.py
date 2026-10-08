"""
The background is subtracted from every well but the reference wells, which are
left out of the result.
"""

from django.test import SimpleTestCase

from plate_calculations.correction import subtract_background


class SubtractBackgroundTest(SimpleTestCase):
    def test_the_background_is_subtracted_from_the_other_wells(self):
        values = {11: 10.0, 12: 4.0, 13: 2.0, 14: 6.0}

        corrected = subtract_background(values, {13, 14}, 4.0)

        self.assertEqual({11: 6.0, 12: 0.0}, corrected)

    def test_reference_wells_without_a_value_are_no_problem(self):
        corrected = subtract_background({11: 10.0}, {13, 14}, 2.0)

        self.assertEqual({11: 8.0}, corrected)
