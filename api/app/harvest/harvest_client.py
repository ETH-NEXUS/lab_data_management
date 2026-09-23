"""
A small client for the Harvest API v2 (https://help.getharvest.com/api-v2/).
"""

import logging

import requests

logger = logging.getLogger(__name__)

# Harvest answers within a second; without a limit a request to a Harvest that
# does not answer would block a server worker forever
TIMEOUT_SECONDS = 10


class HarvestClient:
    def __init__(self, access_token, account_id):
        self.access_token = access_token
        self.account_id = account_id
        self.base_url = "https://api.harvestapp.com/v2/"
        self.headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Harvest-Account-Id": f"{self.account_id}",
            "User-Agent": "Harvest API Example",
            "Content-Type": "application/json",
        }

    def get(self, endpoint, params=None):
        """
        The answer of Harvest as a dict, e.g. get("projects/7") ->
        {"id": 7, "name": "Screening 2026", "notes": "...", ...}

        Raises requests.RequestException when Harvest cannot be reached or
        answers with an error.
        """
        response = requests.get(
            f"{self.base_url}{endpoint}",
            headers=self.headers,
            params=params,
            timeout=TIMEOUT_SECONDS,
        )
        if not response.ok:
            # The answer says what is wrong, e.g. that the token is not valid
            logger.error(
                "Harvest answered %s for %s: %s",
                response.status_code,
                endpoint,
                response.text,
            )
        response.raise_for_status()
        return response.json()
