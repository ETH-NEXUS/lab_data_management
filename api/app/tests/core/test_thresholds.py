"""
Tests for the shared threshold check.

The cases with a zero are the ones that matter most in practice: the Echo
writes `Current Fluid Volume = 0` and `% DMSO = 0` for a failed transfer, and
for years those wells were the only ones that never showed up as problematic.
"""

from django.test import SimpleTestCase

from core.utils.wells.threshold_checks import is_below_threshold, threshold_reasons

THRESHOLD_AMOUNT = 2.5  # microliter
THRESHOLD_DMSO = 80  # percent


class IsBelowThresholdTest(SimpleTestCase):
    def check(self, current_amount, current_dmso):
        return is_below_threshold(
            current_amount, current_dmso, THRESHOLD_AMOUNT, THRESHOLD_DMSO
        )

    def test_a_full_well_is_not_a_problem(self):
        self.assertFalse(self.check(9.5, 95))

    def test_a_low_volume_is_a_problem(self):
        self.assertTrue(self.check(2.0, 95))

    def test_a_low_dmso_is_a_problem(self):
        self.assertTrue(self.check(9.5, 70))

    def test_a_volume_of_zero_is_a_problem(self):
        self.assertTrue(self.check(0, 95))

    def test_a_dmso_of_zero_is_a_problem(self):
        self.assertTrue(self.check(9.5, 0))

    def test_a_failed_transfer_reporting_zeros_is_a_problem(self):
        self.assertTrue(self.check(0, 0))

    def test_a_value_exactly_on_the_volume_threshold_is_not_a_problem(self):
        self.assertFalse(self.check(THRESHOLD_AMOUNT, 95))

    def test_a_value_exactly_on_the_dmso_threshold_is_not_a_problem(self):
        self.assertFalse(self.check(9.5, THRESHOLD_DMSO))

    def test_an_unreported_volume_is_not_a_problem(self):
        self.assertFalse(self.check(None, 95))

    def test_an_unreported_dmso_is_not_a_problem(self):
        self.assertFalse(self.check(9.5, None))

    def test_nothing_reported_at_all_is_not_a_problem(self):
        self.assertFalse(self.check(None, None))

    def test_a_low_volume_counts_even_without_a_dmso_value(self):
        self.assertTrue(self.check(2.0, None))

    def test_a_low_dmso_counts_even_without_a_volume_value(self):
        self.assertTrue(self.check(None, 70))


class ThresholdReasonsTest(SimpleTestCase):
    def reasons(self, current_amount, current_dmso):
        return threshold_reasons(
            current_amount, current_dmso, THRESHOLD_AMOUNT, THRESHOLD_DMSO
        )

    def test_a_full_well_has_no_reason(self):
        self.assertEqual([], self.reasons(9.5, 95))

    def test_a_low_volume_names_the_volume(self):
        self.assertEqual(["volume"], self.reasons(2.0, 95))

    def test_a_low_dmso_names_the_dmso(self):
        self.assertEqual(["dmso"], self.reasons(9.5, 70))

    def test_a_failed_transfer_names_both(self):
        self.assertEqual(["volume", "dmso"], self.reasons(0, 0))

    def test_unreported_values_are_no_reason(self):
        self.assertEqual([], self.reasons(None, None))
