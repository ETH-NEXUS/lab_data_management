"""
The Harvest client. requests.get is replaced by a mock, so no test calls Harvest.
"""

from unittest import mock

import requests
from django.test import SimpleTestCase

from harvest.harvest_client import TIMEOUT_SECONDS, USER_AGENT, HarvestClient


def harvest_response(status_code, text):
    response = requests.Response()
    response.status_code = status_code
    response._content = text.encode()
    return response


class HarvestClientTest(SimpleTestCase):
    def setUp(self):
        self.client = HarvestClient("token", "123")

    @mock.patch("harvest.harvest_client.requests.get")
    def test_the_answer_is_returned_as_a_dict(self, get):
        get.return_value = harvest_response(200, '{"id": 7, "name": "Screening"}')

        answer = self.client.get("projects/7")

        self.assertEqual({"id": 7, "name": "Screening"}, answer)
        get.assert_called_once_with(
            "https://api.harvestapp.com/v2/projects/7",
            headers=self.client.headers,
            params=None,
            timeout=TIMEOUT_SECONDS,
        )

    def test_the_client_names_this_app_to_harvest(self):
        self.assertEqual(USER_AGENT, self.client.headers["User-Agent"])
        self.assertIn("github.com/ETH-NEXUS/lab_data_management", USER_AGENT)

    @mock.patch("harvest.harvest_client.requests.get")
    def test_an_error_of_harvest_is_logged_and_raised(self, get):
        get.return_value = harvest_response(401, '{"error": "invalid_token"}')

        with self.assertLogs("harvest.harvest_client", "ERROR") as logs:
            with self.assertRaises(requests.HTTPError):
                self.client.get("projects")

        self.assertIn("Harvest answered 401 for projects", logs.output[0])
        self.assertIn("invalid_token", logs.output[0])
