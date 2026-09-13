"""
Tests for the endpoint behind the problematic plates card on the messages page.
"""

from datetime import timedelta

from django.contrib.auth import get_user_model
from django.db import connection
from django.test.utils import CaptureQueriesContext
from django.urls import reverse
from django.utils import timezone
from rest_framework import status
from rest_framework.test import APITestCase

from compoundlib.models import CompoundLibrary
from core.models import Plate, PlateDimension, Threshold, Well, WellWithdrawal


class RedFlagViewTest(APITestCase):
    fixtures = ("well_types",)

    def setUp(self):
        Threshold.objects.create(amount=2.5, dmso=80)
        self.dimension = PlateDimension.objects.create(name="dim_4x2", cols=4, rows=2)
        self.library = CompoundLibrary.objects.create(name="Test library")
        self.plate = Plate.objects.create(
            barcode="LIB_PLATE",
            dimension=self.dimension,
            library=self.library,
            status="empty_wells",
        )
        self.user = get_user_model().objects.create(username="tester")
        self.client.force_authenticate(user=self.user)
        self.url = reverse("redflag")

    def add_marked_well(self, position, reports):
        """
        Adds a well that the recalculation has marked, plus what it reported.
        Reports example: [(9.5, 95), (2.0, 95)] - oldest first.
        """
        well = Well.objects.create(plate=self.plate, position=position, status="empty")
        now = timezone.now()
        for index, (current_amount, current_dmso) in enumerate(reports):
            WellWithdrawal.objects.create(
                well=well,
                target_well=None,
                amount=30,
                current_amount=current_amount,
                current_dmso=current_dmso,
                created_at=now - timedelta(days=len(reports) - index),
            )
        return well

    def wells_of_the_plate(self):
        response = self.client.get(self.url)
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        return response.data["Test library"]["LIB_PLATE"]

    def test_anonymous_users_are_turned_away(self):
        self.client.force_authenticate(user=None)
        self.assertEqual(
            status.HTTP_403_FORBIDDEN, self.client.get(self.url).status_code
        )

    def test_a_low_volume_is_reported_with_its_value(self):
        self.add_marked_well(0, [(1.31, 94.5)])
        self.assertEqual(
            [
                {
                    "position": "A1",
                    "current_amount": 1.31,
                    "current_dmso": 94.5,
                    "reasons": ["volume"],
                }
            ],
            self.wells_of_the_plate(),
        )

    def test_a_failed_transfer_reports_both_reasons(self):
        self.add_marked_well(1, [(0, 0)])
        self.assertEqual(["volume", "dmso"], self.wells_of_the_plate()[0]["reasons"])

    def test_a_low_dmso_is_reported_as_such(self):
        self.add_marked_well(2, [(9.5, 70)])
        self.assertEqual(["dmso"], self.wells_of_the_plate()[0]["reasons"])

    def test_the_newest_report_is_the_one_shown(self):
        self.add_marked_well(3, [(9.5, 95), (2.0, 95)])
        entry = self.wells_of_the_plate()[0]
        self.assertEqual(2.0, entry["current_amount"])

    def test_a_well_without_any_report_shows_no_values(self):
        self.add_marked_well(4, [])
        self.assertEqual(
            {
                "position": "B1",
                "current_amount": None,
                "current_dmso": None,
                "reasons": [],
            },
            self.wells_of_the_plate()[0],
        )

    def test_a_flagged_plate_without_marked_wells_stays_in_the_list(self):
        self.assertEqual([], self.wells_of_the_plate())

    def test_plates_that_are_not_flagged_are_left_out(self):
        Plate.objects.filter(id=self.plate.id).update(status=None)
        response = self.client.get(self.url)
        self.assertEqual({}, response.data)

    def test_more_plates_do_not_cost_more_queries(self):
        """The wells of all plates are read in one query, not one per plate."""
        self.add_marked_well(0, [(0, 0)])
        with CaptureQueriesContext(connection) as one_plate:
            self.client.get(self.url)

        for index in range(2):
            plate = Plate.objects.create(
                barcode=f"LIB_PLATE_{index}",
                dimension=self.dimension,
                library=self.library,
                status="empty_wells",
            )
            well = Well.objects.create(plate=plate, position=0, status="empty")
            WellWithdrawal.objects.create(
                well=well, target_well=None, amount=30, current_amount=0, current_dmso=0
            )

        with CaptureQueriesContext(connection) as three_plates:
            self.client.get(self.url)

        self.assertEqual(
            len(one_plate.captured_queries), len(three_plates.captured_queries)
        )
