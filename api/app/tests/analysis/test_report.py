"""
The analysis settings are checked before the report runs: the reports run some of
them as R code.
"""

from django.core.management.base import CommandError
from django.test import SimpleTestCase

from analysis.report import DEFAULT_SETTINGS, check_settings


class CheckSettingsTest(SimpleTestCase):
    def test_settings_that_were_not_chosen_get_the_default(self):
        self.assertEqual(DEFAULT_SETTINGS, check_settings("single", {"act_cut": ""}))

    def test_numbers_from_the_form_become_numbers(self):
        checked = check_settings("single", {"fdr_cut": "0.05", "select_cut": "1"})

        self.assertEqual((0.05, 1.0), (checked["fdr_cut"], checked["select_cut"]))

    def test_a_formula_is_accepted(self):
        checked = check_settings("single", {"act_cut": "log2(1.2) * 3"})

        self.assertEqual("log2(1.2) * 3", checked["act_cut"])

    def test_other_r_code_is_refused(self):
        with self.assertRaisesMessage(CommandError, "act_cut must be a number"):
            check_settings("single", {"act_cut": 'system("ls")'})

    def test_an_unknown_setting_is_refused(self):
        with self.assertRaisesMessage(CommandError, "Unknown analysis setting"):
            check_settings("single", {"path_output": "/etc/"})

    def test_settings_that_are_not_names_with_values_are_refused(self):
        with self.assertRaisesMessage(CommandError, "must be names with values"):
            check_settings("single", "fdr_cut=0.05")

    def test_an_unknown_analysis_type_is_refused(self):
        with self.assertRaisesMessage(CommandError, "Unknown analysis type"):
            check_settings("doseresponse", {})

    def test_a_selectivity_analysis_needs_two_different_conditions(self):
        with self.assertRaisesMessage(CommandError, "two different conditions"):
            check_settings(
                "selectivity", {"condi_yes": "irradiated", "condi_no": "irradiated"}
            )

    def test_a_selectivity_analysis_needs_both_conditions(self):
        with self.assertRaisesMessage(CommandError, "needs both conditions"):
            check_settings("selectivity", {"condi_yes": "irradiated"})
