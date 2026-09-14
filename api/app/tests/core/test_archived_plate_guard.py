"""
Tests that archived library plates and their wells cannot be changed through the API,
while active plates and experiment plates keep working as before.
"""

from unittest.mock import patch

from django.contrib.auth import get_user_model
from rest_framework import status
from rest_framework.test import APITestCase

from compoundlib.models import CompoundLibrary
from core.models import (
    Experiment,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    Well,
    WellDetail,
    WellType,
)

REFUSAL = "Plate ARCHIVED is archived and can no longer be changed."


class ArchivedPlateGuardTest(APITestCase):
    fixtures = ("well_types",)

    def setUp(self):
        dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        library = CompoundLibrary.objects.create(name="Test library")
        self.archived_plate = Plate.objects.create(
            barcode="ARCHIVED", dimension=dimension, library=library, archived=True
        )
        self.active_plate = Plate.objects.create(
            barcode="ACTIVE", dimension=dimension, library=library
        )
        self.archived_well = Well.objects.create(plate=self.archived_plate, position=0)
        self.active_well = Well.objects.create(plate=self.active_plate, position=0)
        user = get_user_model().objects.create(username="tester")
        self.client.force_authenticate(user=user)

    def assert_refused(self, response):
        self.assertEqual(status.HTTP_400_BAD_REQUEST, response.status_code)
        self.assertEqual(REFUSAL, response.data["detail"])

    # Plates

    def test_an_archived_plate_cannot_be_changed(self):
        response = self.client.patch(
            f"/api/plates/{self.archived_plate.id}/", {"barcode": "NEW"}, format="json"
        )
        self.assert_refused(response)
        self.archived_plate.refresh_from_db()
        self.assertEqual("ARCHIVED", self.archived_plate.barcode)

    def test_an_archived_plate_cannot_be_deleted(self):
        response = self.client.delete(f"/api/plates/{self.archived_plate.id}/")
        self.assert_refused(response)
        self.assertTrue(Plate.objects.filter(id=self.archived_plate.id).exists())

    def test_no_template_can_be_applied_to_an_archived_plate(self):
        response = self.client.post(
            f"/api/plates/{self.archived_plate.id}/apply_template/",
            {"template": self.active_plate.id},
            format="json",
        )
        self.assert_refused(response)

    def test_the_archived_flag_itself_can_still_be_changed(self):
        response = self.client.post(
            f"/api/plates/{self.archived_plate.id}/archive/",
            {"archived": False},
            format="json",
        )
        self.assertEqual(status.HTTP_200_OK, response.status_code)

    # Wells

    def test_no_well_can_be_created_on_an_archived_plate(self):
        response = self.client.post(
            "/api/wells/",
            {"plate": self.archived_plate.id, "position": 1},
            format="json",
        )
        self.assert_refused(response)
        self.assertEqual(1, Well.objects.filter(plate=self.archived_plate).count())

    def test_a_well_can_still_be_created_on_an_active_plate(self):
        response = self.client.post(
            "/api/wells/", {"plate": self.active_plate.id, "position": 1}, format="json"
        )
        self.assertEqual(status.HTTP_201_CREATED, response.status_code)

    def test_a_well_of_an_archived_plate_cannot_be_changed(self):
        response = self.client.patch(
            f"/api/wells/{self.archived_well.id}/", {"status": "changed"}, format="json"
        )
        self.assert_refused(response)
        self.archived_well.refresh_from_db()
        self.assertIsNone(self.archived_well.status)

    def test_a_well_cannot_be_moved_to_an_archived_plate(self):
        response = self.client.patch(
            f"/api/wells/{self.active_well.id}/",
            {"plate": self.archived_plate.id, "position": 5},
            format="json",
        )
        self.assert_refused(response)
        self.active_well.refresh_from_db()
        self.assertEqual(self.active_plate.id, self.active_well.plate_id)

    def test_a_well_of_an_archived_plate_cannot_be_deleted(self):
        response = self.client.delete(f"/api/wells/{self.archived_well.id}/")
        self.assert_refused(response)
        self.assertTrue(Well.objects.filter(id=self.archived_well.id).exists())

    def test_a_well_of_an_archived_plate_cannot_be_marked_invalid(self):
        response = self.client.get(
            f"/api/wells/{self.archived_well.id}/mark_as_invalid/"
        )
        self.assert_refused(response)
        self.archived_well.refresh_from_db()
        self.assertFalse(self.archived_well.is_invalid)

    def test_a_well_of_an_active_plate_can_still_be_marked_invalid(self):
        response = self.client.get(f"/api/wells/{self.active_well.id}/mark_as_invalid/")
        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.active_well.refresh_from_db()
        self.assertTrue(self.active_well.is_invalid)

    def test_a_well_of_an_archived_experiment_plate_can_still_be_changed(self):
        """Only library plates become read-only; old experiments keep working."""
        project = Project.objects.create(name="Test project")
        experiment = Experiment.objects.create(name="Test experiment", project=project)
        plate = Plate.objects.create(
            barcode="OLD_EXPERIMENT",
            dimension=self.archived_plate.dimension,
            experiment=experiment,
            archived=True,
        )
        well = Well.objects.create(plate=plate, position=0)
        response = self.client.get(f"/api/wells/{well.id}/mark_as_invalid/")
        self.assertEqual(status.HTTP_200_OK, response.status_code)

    # The response reads materialized views, which the test database does not refresh.
    @patch("core.views.plates.PlateSerializer")
    @patch.object(WellDetail, "refresh")
    @patch.object(PlateDetail, "refresh")
    def test_a_template_for_all_plates_skips_archived_plates(
        self, plate_refresh, well_refresh, plate_serializer
    ):
        plate_serializer.return_value.data = {}
        template = Plate.objects.create(
            barcode="TEMPLATE", dimension=self.archived_plate.dimension
        )
        Well.objects.create(plate=template, position=0, type=WellType.by_name("P"))

        response = self.client.post(
            f"/api/plates/{self.active_plate.id}/apply_template/",
            {"template": template.id, "apply_to_all_experiment_plates": True},
            format="json",
        )

        self.assertEqual(status.HTTP_200_OK, response.status_code)
        self.active_well.refresh_from_db()
        self.assertEqual("P", self.active_well.type.name)
        self.archived_well.refresh_from_db()
        self.assertNotEqual("P", self.archived_well.type.name)
        self.assertEqual(1, Well.objects.filter(plate=self.archived_plate).count())
