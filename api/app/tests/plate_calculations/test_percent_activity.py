"""
%Activity = 100 * (value - median(P)) / (median(N) - median(P)), from the raw
readout and per time point; the negative control wells are left empty.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase
from django.urls import reverse

from plate_calculations.percent_activity import percent_activity
from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDimension,
    Project,
    Well,
    WellType,
)

FIRST_READ = datetime(2026, 9, 30, 15, 48)
SECOND_READ = datetime(2026, 9, 30, 17, 0)


class PercentActivityTest(SimpleTestCase):
    def test_n_is_100_and_p_is_0_percent(self):
        # Medians: N 200, P 10
        values = {11: 100.0, 12: 300.0, 21: 0.0, 22: 20.0, 31: 105.0}

        activity = percent_activity(values, {11, 12}, 200.0, 10.0)

        self.assertEqual([21, 22, 31], sorted(activity))
        self.assertAlmostEqual(-100 * 10 / 190, activity[21])
        self.assertAlmostEqual(100 * 10 / 190, activity[22])
        # Halfway between P and N
        self.assertAlmostEqual(50.0, activity[31])


class PercentActivityViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_1x5", rows=1, cols=5)
        self.plate = Plate.objects.create(
            barcode="RKS_1", dimension=dimension, experiment=experiment
        )
        negative = WellType.objects.create(name="N", description="negative")
        positive = WellType.objects.create(name="P", description="positive")
        compound = WellType.objects.create(name="C", description="compound")
        # Well type and the values of the first and the second read
        layout = [
            (negative, 2900.0, 5800.0),
            (negative, 3000.0, 6000.0),
            (positive, 600.0, 1200.0),
            (positive, 700.0, 1400.0),
            (compound, 1800.0, 3600.0),
        ]
        for position, (well_type, first, second) in enumerate(layout):
            well = Well.objects.create(
                plate=self.plate, position=position, type=well_type
            )
            Measurement.objects.create(
                well=well, label="Lum1", value=first, measured_at=FIRST_READ
            )
            Measurement.objects.create(
                well=well, label="Lum1", value=second, measured_at=SECOND_READ
            )

    def activity(self, data):
        url = reverse("percent_activity_of_plate", args=[self.plate.id])
        return self.client.post(url, data, content_type="application/json")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def test_without_login_nothing_is_saved(self):
        response = self.activity(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        self.assertEqual(403, response.status_code)
        self.assertFalse(Measurement.objects.exclude(label="Lum1").exists())

    def test_the_activity_of_every_read_is_saved(self):
        self.login()

        response = self.activity(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        self.assertEqual({"label": "Lum1_activity_N_P"}, response.json())
        saved = {
            (position, measured_at): value
            for position, measured_at, value in Measurement.objects.filter(
                label="Lum1_activity_N_P"
            ).values_list("well__position", "measured_at", "value")
        }
        # The N wells are left empty; medians N 2950 / 5900, P 650 / 1300
        self.assertEqual(
            {(2, FIRST_READ), (3, FIRST_READ), (4, FIRST_READ)},
            {key for key in saved if key[1] == FIRST_READ},
        )
        self.assertAlmostEqual(100 * 1150 / 2300, saved[(4, FIRST_READ)])
        self.assertAlmostEqual(100 * 2300 / 4600, saved[(4, SECOND_READ)])
        self.assertAlmostEqual(100 * -50 / 2300, saved[(2, FIRST_READ)])

    def test_the_same_well_type_for_both_controls_is_refused(self):
        self.login()

        response = self.activity(
            {"label": "Lum1", "negative_type": "N", "positive_type": "N"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("must be different", str(response.json()))

    def test_a_missing_control_is_refused(self):
        self.login()

        response = self.activity(
            {"label": "Lum1", "negative_type": "N1", "positive_type": "P"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn('no "N1" wells', response.json()[0])

    def test_controls_with_the_same_median_are_refused(self):
        self.login()
        Measurement.objects.filter(label="Lum1").update(value=500.0)

        response = self.activity(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("are the same", response.json()[0])
