"""
The normalization of the R report (Michael's analysis), for every well:
%Inhibition = (log10(x) - median(log10 N)) / (median(log10 P) - median(log10 N))
and %Activity = 1 - %Inhibition, as fractions and per time point.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase
from django.urls import reverse

from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDimension,
    Project,
    Well,
    WellType,
)
from plate_calculations.normalization import inhibition

FIRST_READ = datetime(2026, 9, 30, 15, 48)
SECOND_READ = datetime(2026, 9, 30, 17, 0)


class InhibitionTest(SimpleTestCase):
    def test_n_is_0_and_p_is_1(self):
        logs = {11: 3.0, 21: 1.0, 31: 2.0}

        self.assertEqual({11: 0.0, 21: 1.0, 31: 0.5}, inhibition(logs, 3.0, 1.0))


class NormalizationViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_1x6", rows=1, cols=6)
        self.plate = Plate.objects.create(
            barcode="RKS_1", dimension=dimension, experiment=experiment
        )
        negative = WellType.objects.create(name="N", description="negative")
        positive = WellType.objects.create(name="P", description="positive")
        compound = WellType.objects.create(name="C", description="compound")
        # Well type and the values of the first and the second read; log10 of
        # N 1000 = 3, P 10 = 1, so the compound with 100 (log10 2) is halfway
        layout = [
            (negative, 1000.0, 10000.0),
            (negative, 1000.0, 10000.0),
            (positive, 10.0, 100.0),
            (positive, 10.0, 100.0),
            (compound, 100.0, 1000.0),
            (compound, 0.0, 50.0),
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

    def normalize(self, data):
        url = reverse("normalize_plate_measurement", args=[self.plate.id])
        return self.client.post(url, data, content_type="application/json")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def saved(self, label):
        """{(position, measured_at): value} of one measurement."""
        return {
            (position, measured_at): value
            for position, measured_at, value in Measurement.objects.filter(
                label=label
            ).values_list("well__position", "measured_at", "value")
        }

    def test_without_login_nothing_is_saved(self):
        response = self.normalize(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        self.assertEqual(403, response.status_code)
        self.assertFalse(Measurement.objects.exclude(label="Lum1").exists())

    def test_inhibition_and_activity_of_every_well_and_read(self):
        self.login()

        response = self.normalize(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        # The well with 0 has no log10: it is left empty in every read
        self.assertEqual(
            {
                "label": "Lum1_inhibition_N_P",
                "activity_label": "Lum1_activity_N_P",
                "skipped": 1,
            },
            response.json(),
        )
        inhibition = self.saved("Lum1_inhibition_N_P")
        activity = self.saved("Lum1_activity_N_P")
        for read in [FIRST_READ, SECOND_READ]:
            self.assertEqual(0.0, inhibition[(0, read)])
            self.assertEqual(1.0, inhibition[(2, read)])
            self.assertAlmostEqual(0.5, inhibition[(4, read)])
            self.assertEqual(1.0, activity[(0, read)])
            self.assertAlmostEqual(0.5, activity[(4, read)])
        self.assertEqual(10, len(inhibition))
        self.assertEqual(10, len(activity))

    def test_the_same_well_type_for_both_controls_is_refused(self):
        self.login()

        response = self.normalize(
            {"label": "Lum1", "negative_type": "N", "positive_type": "N"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("must be different", str(response.json()))

    def test_a_missing_control_is_refused(self):
        self.login()

        response = self.normalize(
            {"label": "Lum1", "negative_type": "N1", "positive_type": "P"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn('no "N1" wells', response.json()[0])

    def test_controls_with_the_same_median_are_refused(self):
        self.login()
        Measurement.objects.filter(label="Lum1").update(value=500.0)

        response = self.normalize(
            {"label": "Lum1", "negative_type": "N", "positive_type": "P"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("cannot be normalized", response.json()[0])
        self.assertFalse(Measurement.objects.exclude(label="Lum1").exists())
