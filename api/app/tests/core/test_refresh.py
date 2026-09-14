"""
Tests for the endpoint behind the refresh button in the page header.

Refreshing rewrites the materialized views of the whole database, so it must
only run on a POST from a logged in user that carries a CSRF token. The refresh
itself is replaced by a stand-in here: these tests are about who may start it.
"""

from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient, APITestCase

from core.models import ExperimentDetail, PlateDetail, WellDetail


@patch.object(ExperimentDetail, "refresh")
@patch.object(WellDetail, "refresh")
@patch.object(PlateDetail, "refresh")
class RefreshTest(APITestCase):
    def setUp(self):
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.url = reverse("refresh")

    def test_a_logged_in_user_can_refresh(
        self, plate_refresh, well_refresh, experiment_refresh
    ):
        self.client.force_authenticate(user=self.user)
        response = self.client.post(self.url)
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.assertEqual({"status": "Data refreshed successfully"}, response.data)
        plate_refresh.assert_called_once_with(concurrently=True)
        well_refresh.assert_called_once_with(concurrently=True)
        experiment_refresh.assert_called_once_with(concurrently=True)

    def test_get_is_not_allowed(self, plate_refresh, well_refresh, experiment_refresh):
        self.client.force_authenticate(user=self.user)
        response = self.client.get(self.url)
        self.assertEqual(status.HTTP_405_METHOD_NOT_ALLOWED, response.status_code)
        plate_refresh.assert_not_called()

    def test_anonymous_users_are_turned_away(
        self, plate_refresh, well_refresh, experiment_refresh
    ):
        response = self.client.post(self.url)
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)
        plate_refresh.assert_not_called()

    def test_a_session_without_a_csrf_token_is_turned_away(
        self, plate_refresh, well_refresh, experiment_refresh
    ):
        # A real browser session, where DRF checks the CSRF token.
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        response = client.post(self.url)
        self.assertEqual(status.HTTP_403_FORBIDDEN, response.status_code)
        plate_refresh.assert_not_called()

    def test_a_session_with_a_csrf_token_can_refresh(
        self, plate_refresh, well_refresh, experiment_refresh
    ):
        client = APIClient(enforce_csrf_checks=True)
        client.login(username="tester", password="test-password")
        client.get(reverse("auth-cookie"))
        csrf_token = client.cookies["csrftoken"].value

        response = client.post(self.url, HTTP_X_CSRFTOKEN=csrf_token)
        self.assertEqual(status.HTTP_200_OK, response.status_code)

    def test_a_failure_does_not_reveal_details(
        self, plate_refresh, well_refresh, experiment_refresh
    ):
        plate_refresh.side_effect = RuntimeError("secret database detail")
        self.client.force_authenticate(user=self.user)
        with self.assertLogs("core.views.system", level="ERROR"):
            response = self.client.post(self.url)
        self.assertEqual(status.HTTP_500_INTERNAL_SERVER_ERROR, response.status_code)
        self.assertNotIn("secret database detail", str(response.data))
