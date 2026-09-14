"""
Tests for creating a plate mapping through the API, which copies one plate onto another.
"""

from unittest.mock import patch

from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.test import APITestCase

from compoundlib.models import Compound
from core.models import (
    Plate,
    PlateDetail,
    PlateDimension,
    PlateMapping,
    Well,
    WellCompound,
    WellDetail,
)


# The materialized views are refreshed after every mapping; that is not what is tested here.
@patch.object(WellDetail, "refresh")
@patch.object(PlateDetail, "refresh")
class PlateMappingApiTest(APITestCase):
    fixtures = ("well_types",)

    def setUp(self):
        dimension = PlateDimension.objects.create(name="dim_2x1", rows=1, cols=2)
        self.source_plate = Plate.objects.create(barcode="SOURCE", dimension=dimension)
        self.target_plate = Plate.objects.create(barcode="TARGET", dimension=dimension)
        compound = Compound.objects.create(name="compound_a", structure="C")
        source_well = Well.objects.create(plate=self.source_plate, position=1)
        WellCompound.objects.create(well=source_well, compound=compound, amount=10)
        self.user = get_user_model().objects.create_user(
            username="tester", password="test-password"
        )
        self.client.force_authenticate(user=self.user)

    def test_a_plate_copy_creates_the_mapping_and_copies_the_wells(
        self, plate_refresh, well_refresh
    ):
        response = self.client.post(
            "/api/platemappings/",
            {
                "source_plate": self.source_plate.id,
                "target_plate": self.target_plate.id,
                "amount": 5,
                # The web page sends "undefined" for empty fields; it is stored as null.
                "from_column": "undefined",
            },
            format="json",
        )

        self.assertEqual(status.HTTP_201_CREATED, response.status_code, response.content)
        mapping = PlateMapping.objects.get()
        self.assertIsNone(mapping.from_column)
        self.assertEqual(5, mapping.amount)

        target_well = self.target_plate.wells.get(position=1)
        self.assertEqual("compound_a", target_well.compounds.get().name)
        self.assertEqual(5, target_well.donors.get().amount)
        plate_refresh.assert_called()
        well_refresh.assert_called()
