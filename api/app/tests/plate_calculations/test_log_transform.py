"""
The log10 of a measurement of one plate is saved as a new measurement; values
of 0 or below are left out and counted. Only for logged in users.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase
from django.urls import reverse

from plate_calculations.log_transform import wells_without_log10
from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    Well,
    WellType,
)

FIRST_READ = datetime(2026, 9, 30, 15, 48)
SECOND_READ = datetime(2026, 9, 30, 17, 0)


class WellsWithoutLog10Test(SimpleTestCase):
    def test_a_value_of_0_or_below_in_any_read_leaves_the_well_out(self):
        values = {
            FIRST_READ: {11: 1000.0, 12: 0.0, 13: 5.0},
            SECOND_READ: {11: 10.0, 12: 5.0, 13: -1.0},
        }

        self.assertEqual({12, 13}, wells_without_log10(values))


class Log10ViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_1x3", rows=1, cols=3)
        self.plate = Plate.objects.create(
            barcode="RKS_1", dimension=dimension, experiment=experiment
        )
        well_type = WellType.objects.create(name="R", description="reference")
        # Values of the first and the second read
        self.wells = []
        for position, (first, second) in enumerate(
            [(100.0, 1000.0), (0.0, 10.0), (1.0, -3.0)]
        ):
            well = Well.objects.create(
                plate=self.plate, position=position, type=well_type
            )
            self.wells.append(well)
            Measurement.objects.create(
                well=well, label="Lum1", value=first, measured_at=FIRST_READ
            )
            Measurement.objects.create(
                well=well, label="Lum1", value=second, measured_at=SECOND_READ
            )

    def log10(self, data, plate_id=None):
        url = reverse("log10_plate_measurement", args=[plate_id or self.plate.id])
        return self.client.post(url, data, content_type="application/json")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def test_without_login_nothing_is_saved(self):
        self.assertEqual(403, self.log10({"label": "Lum1"}).status_code)
        self.assertFalse(Measurement.objects.filter(label="Lum1_log10").exists())

    def test_the_log10_of_every_read_is_saved(self):
        self.login()

        response = self.log10({"label": "Lum1"})

        # Well 1 has 0 in the first read, well 2 has -3 in the second one:
        # both are left empty in every read
        self.assertEqual({"label": "Lum1_log10", "skipped": 2}, response.json())
        saved = sorted(
            Measurement.objects.filter(label="Lum1_log10").values_list(
                "well__position", "measured_at", "value"
            )
        )
        self.assertEqual([(0, FIRST_READ, 2.0), (0, SECOND_READ, 3.0)], saved)
        self.assertIn(
            "Lum1_log10", PlateDetail.objects.get(id=self.plate.id).measurement_labels
        )

    def test_a_new_log10_replaces_the_last_one(self):
        self.login()
        self.log10({"label": "Lum1"})

        self.log10({"label": "Lum1"})

        self.assertEqual(2, Measurement.objects.filter(label="Lum1_log10").count())

    def test_a_measurement_without_values_above_0_is_refused(self):
        self.login()
        Measurement.objects.filter(label="Lum1").update(value=0.0)

        response = self.log10({"label": "Lum1"})

        self.assertEqual(400, response.status_code)
        self.assertIn("has a value of 0 or below", response.json()[0])

    def test_an_unknown_measurement_is_refused(self):
        self.login()

        response = self.log10({"label": "Fluo"})

        self.assertEqual(400, response.status_code)
        self.assertIn('no measurement "Fluo"', response.json()[0])
