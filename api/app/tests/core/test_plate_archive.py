"""
Tests for archiving and unarchiving a plate through the API, as the plate page does.
"""

from django.contrib.auth import get_user_model
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient, APITestCase

from core.models import Plate


class PlateArchiveTest(APITestCase):
    def setUp(self):
        self.plate = Plate.objects.create(barcode="ARCHIVE_ME")
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.url = f"/api/plates/{self.plate.id}/archive/"

    def archived_in_database(self):
        self.plate.refresh_from_db()
        return self.plate.archived

    def test_a_logged_in_user_can_archive_a_plate(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.post(self.url, {"archived": True}, format="json")
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertEqual(
            {"id": self.plate.id, "barcode": "ARCHIVE_ME", "archived": True},
            response.data,
        )
        self.assertTrue(self.archived_in_database())

    def test_a_logged_in_user_can_unarchive_a_plate(self):
        Plate.objects.filter(pk=self.plate.pk).update(archived=True)
        self.client.force_authenticate(user=self.user)
        response = self.client.post(self.url, {"archived": False}, format="json")
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertFalse(self.archived_in_database())

    def test_anything_but_true_or_false_is_refused(self):
        self.client.force_authenticate(user=self.user)
        for body in ({"archived": "yes"}, {"archived": 1}, {"archived": None}, {}):
            response = self.client.post(self.url, body, format="json")
            self.assertEqual(status.HTTP_400_BAD_REQUEST, response.status_code, body)
        self.assertFalse(self.archived_in_database())

    def test_anonymous_users_are_turned_away(self):
        response = self.client.post(self.url, {"archived": True}, format="json")
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)
        self.assertFalse(self.archived_in_database())

    def test_get_is_not_allowed(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.get(self.url)
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)

    def test_an_unknown_plate_is_not_found(self):
        self.client.force_authenticate(user=self.user)
        for plate_id in (self.plate.id + 1000, "not-a-number"):
            response = self.client.post(
                f"/api/plates/{plate_id}/archive/", {"archived": True}, format="json"
            )
            self.assertEqual(status.HTTP_404_NOT_FOUND, response.status_code, plate_id)

    def test_a_session_without_a_csrf_token_is_turned_away(self):
        # A real browser session, where DRF checks the CSRF token.
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        response = client.post(self.url, {"archived": True}, format="json")
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)
        self.assertFalse(self.archived_in_database())

    def test_a_session_with_a_csrf_token_can_archive(self):
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        client.get(reverse("auth-cookie"))
        csrf_token = client.cookies["csrftoken"].value

        response = client.post(
            self.url, {"archived": True}, format="json", HTTP_X_CSRFTOKEN=csrf_token
        )
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertTrue(self.archived_in_database())
