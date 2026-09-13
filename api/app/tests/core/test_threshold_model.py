"""
Tests for the one threshold of the application: how it is read, and which
values it accepts.
"""

from django.contrib.auth import get_user_model
from django.core.exceptions import ValidationError
from django.core.management import call_command
from django.test import TestCase
from rest_framework import status
from rest_framework.test import APITestCase

from compoundlib.models import CompoundLibrary
from core.models import Plate, PlateDimension, Threshold, Well, WellWithdrawal


class ThresholdCurrentTest(TestCase):
    def test_the_existing_threshold_is_returned(self):
        existing = Threshold.objects.create(pk=5, amount=3.0, dmso=70)
        self.assertEqual(existing, Threshold.current())
        self.assertEqual(1, Threshold.objects.count())

    def test_a_missing_threshold_is_created_with_the_defaults(self):
        threshold = Threshold.current()
        self.assertEqual(2.5, threshold.amount)
        self.assertEqual(80, threshold.dmso)

    def test_asking_twice_does_not_create_a_second_threshold(self):
        Threshold.current()
        Threshold.current()
        self.assertEqual(1, Threshold.objects.count())


class ThresholdValidationTest(TestCase):
    def assert_rejected(self, **values):
        with self.assertRaises(ValidationError):
            Threshold(**values).full_clean()

    def test_a_negative_volume_is_rejected(self):
        self.assert_rejected(amount=-1, dmso=80)

    def test_a_negative_dmso_is_rejected(self):
        self.assert_rejected(amount=2.5, dmso=-1)

    def test_a_dmso_above_100_percent_is_rejected(self):
        self.assert_rejected(amount=2.5, dmso=101)

    def test_the_limits_themselves_are_accepted(self):
        Threshold(amount=0, dmso=0).full_clean()
        Threshold(amount=2.5, dmso=100).full_clean()


class ThresholdApiValidationTest(APITestCase):
    """The form on the messages page saves through this endpoint."""

    def setUp(self):
        self.threshold = Threshold.objects.create(amount=2.5, dmso=80)
        user = get_user_model().objects.create(username="tester")
        self.client.force_authenticate(user=user)
        self.url = f"/api/thresholds/{self.threshold.id}/"

    def test_a_valid_change_is_saved(self):
        response = self.client.patch(self.url, {"amount": 2.0, "dmso": 75})
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.threshold.refresh_from_db()
        self.assertEqual(2.0, self.threshold.amount)

    def test_a_dmso_above_100_percent_is_refused(self):
        response = self.client.patch(self.url, {"dmso": 101})
        self.assertEqual(status.HTTP_400_BAD_REQUEST, response.status_code)

    def test_a_negative_volume_is_refused(self):
        response = self.client.patch(self.url, {"amount": -1})
        self.assertEqual(status.HTTP_400_BAD_REQUEST, response.status_code)


class RecalculationWithoutThresholdRowTest(TestCase):
    """Without a threshold row the defaults apply, instead of marking nothing."""

    fixtures = ("well_types",)

    def test_an_empty_well_is_still_marked(self):
        dimension = PlateDimension.objects.create(name="dim_2x1", cols=2, rows=1)
        library = CompoundLibrary.objects.create(name="Test library")
        plate = Plate.objects.create(
            barcode="LIB_PLATE", dimension=dimension, library=library
        )
        well = Well.objects.create(plate=plate, position=0)
        WellWithdrawal.objects.create(
            well=well, target_well=None, amount=30, current_amount=0, current_dmso=95
        )

        call_command("find_problems", "mark_empty_wells", verbosity=0)

        well.refresh_from_db()
        self.assertEqual("empty", well.status)
