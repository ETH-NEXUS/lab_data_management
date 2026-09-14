"""
Tests for the check of the request data when a template is applied to a plate.
"""

from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.test import APITestCase

from core.models import Plate, PlateDimension


class ApplyTemplateValidationTest(APITestCase):
    def setUp(self):
        dimension = PlateDimension.objects.create(name="dim_2x2", rows=2, cols=2)
        self.plate = Plate.objects.create(barcode="NEEDS_TEMPLATE", dimension=dimension)
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.client.force_authenticate(user=self.user)

    def test_a_missing_template_is_refused_with_400(self):
        # This used to answer 500: the view raised Django's ValidationError,
        # which the REST framework does not turn into a response.
        response = self.client.post(
            f"/api/plates/{self.plate.id}/apply_template/", {}, format="json"
        )
        self.assertEqual(status.HTTP_400_BAD_REQUEST, response.status_code)
        self.assertEqual({"template": ["This field is required."]}, response.json())
        self.assertEqual(0, self.plate.wells.count())
