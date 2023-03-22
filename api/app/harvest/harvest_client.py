import requests


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
        response = requests.get(
            f"{self.base_url}{endpoint}", headers=self.headers, params=params
        )
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            print("Error:", e)
            print("Response headers:", response.headers)
            print("Response body:", response.text)
            raise
        return response.json()
