"""
Deleting a calculated measurement of one plate from the plate page: all its
values go, an imported measurement is refused. Only for logged in users.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from compoundlib.models import CompoundLibrary
from core.models import (
    Experiment,
    Measurement,
    MeasurementAssignment,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    Well,
    WellType,
)

FIRST_READ = datetime(2026, 9, 30, 15, 48)
SECOND_READ = datetime(2026, 9, 30, 17, 0)


class DeletionViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_1x2", rows=1, cols=2)
        self.plate = Plate.objects.create(
            barcode="RKS_1", dimension=dimension, experiment=experiment
        )
        well_type = WellType.objects.create(name="N", description="negative")
        self.wells = []
        for position in range(2):
            well = Well.objects.create(
                plate=self.plate, position=position, type=well_type
            )
            self.wells.append(well)
            for read in [FIRST_READ, SECOND_READ]:
                Measurement.objects.create(
                    well=well, label="Lum1", value=100.0, measured_at=read
                )
                Measurement.objects.create(
                    well=well, label="Lum1_log10", value=2.0, measured_at=read
                )

    def delete(self, data):
        url = reverse("delete_plate_calculation", args=[self.plate.id])
        return self.client.post(url, data, content_type="application/json")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def labels(self):
        return set(
            Measurement.objects.filter(well__plate=self.plate).values_list(
                "label", flat=True
            )
        )

    def test_without_login_nothing_is_deleted(self):
        self.assertEqual(403, self.delete({"label": "Lum1_log10"}).status_code)
        self.assertEqual({"Lum1", "Lum1_log10"}, self.labels())

    def test_every_value_of_the_measurement_is_deleted(self):
        self.login()

        response = self.delete({"label": "Lum1_log10"})

        self.assertEqual(200, response.status_code)
        # 2 wells, 2 reads
        self.assertEqual({"label": "Lum1_log10", "deleted": 4}, response.json())
        self.assertEqual({"Lum1"}, self.labels())
        # The plate page reads the refreshed view
        details = PlateDetail.objects.get(pk=self.plate.id)
        self.assertEqual(["Lum1"], details.measurement_labels)

    def test_an_imported_measurement_is_not_deleted(self):
        self.login()
        assignment = MeasurementAssignment.objects.create(
            plate=self.plate, filename="RKS_1.asc"
        )
        Measurement.objects.filter(label="Lum1").update(
            measurement_assignment=assignment
        )

        response = self.delete({"label": "Lum1"})

        self.assertEqual(400, response.status_code)
        self.assertIn("imported from a file", response.json()[0])
        self.assertEqual({"Lum1", "Lum1_log10"}, self.labels())

    def test_an_unknown_measurement_is_refused(self):
        self.login()

        response = self.delete({"label": "Fluo"})

        self.assertEqual(400, response.status_code)
        self.assertIn('no measurement "Fluo"', response.json()[0])

    def test_an_archived_library_plate_is_not_changed(self):
        self.login()
        # A plate belongs to a library or to an experiment, not to both
        self.plate.experiment = None
        self.plate.library = CompoundLibrary.objects.create(name="Library")
        self.plate.archived = True
        self.plate.save()

        response = self.delete({"label": "Lum1_log10"})

        self.assertEqual(400, response.status_code)
        self.assertIn("is archived", response.json()["detail"])
        self.assertEqual({"Lum1", "Lum1_log10"}, self.labels())
