"""
Tests for the endpoint behind the "Recalculate" button on the messages page.

The recalculation writes to the database, so it must only run on a POST from a
logged in user that carries a CSRF token.
"""

from django.contrib.auth import get_user_model
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient, APITestCase


class RecalculateStatusTest(APITestCase):
    def setUp(self):
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.url = reverse("recalculate_status")

    def test_a_logged_in_user_can_recalculate(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.post(self.url)
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertEqual({"status": "ok"}, response.data)

    def test_get_is_not_allowed(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.get(self.url)
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)

    def test_anonymous_users_are_turned_away(self):
        response = self.client.post(self.url)
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)

    def test_a_session_without_a_csrf_token_is_turned_away(self):
        # A real browser session, where DRF checks the CSRF token.
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        response = client.post(self.url)
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)

    def test_a_session_with_a_csrf_token_can_recalculate(self):
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        client.get(reverse("auth-cookie"))
        csrf_token = client.cookies["csrftoken"].value

        response = client.post(self.url, HTTP_X_CSRFTOKEN=csrf_token)
        self.assertEqual(status.HTTP_200_OK, response.status_code)
