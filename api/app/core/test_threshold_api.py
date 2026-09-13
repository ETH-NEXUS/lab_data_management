"""
Tests for what the threshold API allows. The messages page reads the threshold
and changes it; nothing else should be possible through the API.
"""

from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.test import APIClient, APITestCase

from core.models import Threshold


class ThresholdApiAccessTest(APITestCase):
    def setUp(self):
        self.threshold = Threshold.objects.create(amount=2.5, dmso=80)
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.client.force_authenticate(user=self.user)
        self.list_url = "/api/thresholds/"
        self.detail_url = f"/api/thresholds/{self.threshold.id}/"

    def test_anonymous_users_are_turned_away(self):
        self.client.force_authenticate(user=None)
        response = self.client.get(self.list_url)
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)

    def test_the_list_keeps_the_shape_the_page_reads(self):
        response = self.client.get(self.list_url)
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertEqual(self.threshold.id, response.data["results"][0]["id"])

    def test_a_single_threshold_can_be_read(self):
        response = self.client.get(self.detail_url)
        self.assertEqual(status.HTTP_200_OK, response.status_code)

    def test_a_logged_in_user_can_change_the_threshold(self):
        response = self.client.patch(self.detail_url, {"amount": 2.0})
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.threshold.refresh_from_db()
        self.assertEqual(2.0, self.threshold.amount)

    def test_a_second_threshold_cannot_be_created(self):
        response = self.client.post(self.list_url, {"amount": 1.0, "dmso": 50})
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)
        self.assertEqual(1, Threshold.objects.count())

    def test_the_threshold_cannot_be_deleted(self):
        response = self.client.delete(self.detail_url)
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)
        self.assertTrue(Threshold.objects.filter(id=self.threshold.id).exists())

    def test_the_threshold_cannot_be_replaced_with_put(self):
        response = self.client.put(self.detail_url, {"amount": 1.0, "dmso": 50})
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)

    def test_a_change_without_a_csrf_token_is_turned_away(self):
        # A real browser session, where DRF checks the CSRF token.
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        response = client.patch(self.detail_url, {"amount": 2.0})
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)
