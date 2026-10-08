"""
The background correction of one plate is started from the plate page and saved
as a new measurement of the plate, which the plate page reads from the
materialized views. Only for logged in users.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from core.models import (
    Experiment,
    Measurement,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    Well,
    WellDetail,
    WellType,
)

FIRST_READ = datetime(2025, 5, 16, 10, 0)
SECOND_READ = datetime(2025, 5, 16, 12, 0)


class BackgroundCorrectionViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_2x3", rows=2, cols=3)
        self.plate = Plate.objects.create(
            barcode="SP_1", dimension=dimension, experiment=experiment
        )
        reference = WellType.objects.create(name="Nref", description="reference")
        positive = WellType.objects.create(name="P", description="positive")
        compound = WellType.objects.create(name="C", description="compound")

        # Well type and the values of the first and the second read
        layout = [
            (reference, 2.0, 20.0),
            (reference, 4.0, 40.0),
            (reference, 9.0, 90.0),
            (positive, 1.0, 10.0),
            (compound, 10.0, 100.0),
            (compound, 7.0, 70.0),
        ]
        self.wells = []
        for position, (well_type, first, second) in enumerate(layout):
            well = Well.objects.create(
                plate=self.plate, position=position, type=well_type
            )
            self.wells.append(well)
            Measurement.objects.create(
                well=well, label="Lum", value=first, measured_at=FIRST_READ
            )
            Measurement.objects.create(
                well=well, label="Lum", value=second, measured_at=SECOND_READ
            )

    def correct(self, settings):
        url = reverse("correct_plate_background", args=[self.plate.id])
        return self.client.post(url, settings, content_type="application/json")

    def login(self):
        self.client.force_login(User.objects.create_user("tester"))

    def test_without_login_nothing_is_saved(self):
        response = self.correct(
            {"label": "Lum", "reference_type": "Nref", "method": "median"}
        )

        self.assertEqual(403, response.status_code)
        self.assertFalse(Measurement.objects.exclude(label="Lum").exists())

    def test_median_is_subtracted_per_time_point_and_shown_on_the_plate(self):
        self.login()

        response = self.correct(
            {"label": "Lum", "reference_type": "Nref", "method": "median"}
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual({"label": "Lum_bc_Nref_median"}, response.json())
        # Median of the reference wells: 4 at the first read, 40 at the second
        well_details = WellDetail.objects.filter(plate_id=self.plate.id)
        corrected = {
            well.position: well.measurements.get("Lum_bc_Nref_median")
            for well in well_details
        }
        self.assertEqual(
            {
                0: None,
                1: None,
                2: None,
                3: [-3.0, -30.0],
                4: [6.0, 60.0],
                5: [3.0, 30.0],
            },
            corrected,
        )
        plate_detail = PlateDetail.objects.get(id=self.plate.id)
        self.assertIn("Lum_bc_Nref_median", plate_detail.measurement_labels)
        # The reference wells have no corrected value, so no stats either
        self.assertEqual(
            ["C", "P"], sorted(plate_detail.stats["Lum_bc_Nref_median"].keys())
        )

    def test_mean_gets_its_own_measurement(self):
        self.login()

        self.correct({"label": "Lum", "reference_type": "Nref", "method": "median"})
        response = self.correct(
            {"label": "Lum", "reference_type": "Nref", "method": "mean"}
        )

        self.assertEqual({"label": "Lum_bc_Nref_mean"}, response.json())
        # Mean of the reference wells at the first read: (2 + 4 + 9) / 3 = 5
        compound = Measurement.objects.get(
            well=self.wells[4], label="Lum_bc_Nref_mean", measured_at=FIRST_READ
        )
        self.assertEqual(5.0, compound.value)
        self.assertTrue(Measurement.objects.filter(label="Lum_bc_Nref_median").exists())

    def test_a_new_correction_replaces_the_last_one(self):
        self.login()
        self.correct({"label": "Lum", "reference_type": "Nref", "method": "median"})
        Measurement.objects.filter(
            well=self.wells[1], label="Lum", measured_at=FIRST_READ
        ).update(value=5.0)

        self.correct({"label": "Lum", "reference_type": "Nref", "method": "median"})

        # Median of 2, 5 and 9 is now 5
        compound = Measurement.objects.get(
            well=self.wells[4], label="Lum_bc_Nref_median", measured_at=FIRST_READ
        )
        self.assertEqual(5.0, compound.value)
        # 3 wells x 2 reads, no copies of the first run
        self.assertEqual(
            6, Measurement.objects.filter(label="Lum_bc_Nref_median").count()
        )

    def test_reference_type_that_is_not_on_the_plate_is_refused(self):
        self.login()

        response = self.correct(
            {"label": "Lum", "reference_type": "N1", "method": "median"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn('no "N1" wells', response.json()[0])
        self.assertFalse(Measurement.objects.exclude(label="Lum").exists())

    def test_label_ending_with_a_space_is_corrected(self):
        self.login()
        Measurement.objects.filter(label="Lum").update(label="Lum ")

        response = self.correct(
            {"label": "Lum ", "reference_type": "Nref", "method": "median"}
        )

        self.assertEqual({"label": "Lum _bc_Nref_median"}, response.json())

    def test_a_result_without_values_is_refused_and_keeps_the_last_one(self):
        self.login()
        self.correct({"label": "Lum", "reference_type": "Nref", "method": "median"})
        # Only the reference wells still have a value, so nothing is left to correct
        Measurement.objects.filter(label="Lum").exclude(
            well__type__name="Nref"
        ).delete()

        response = self.correct(
            {"label": "Lum", "reference_type": "Nref", "method": "median"}
        )

        self.assertEqual(400, response.status_code)
        self.assertIn("gets a value", response.json()[0])
        self.assertEqual(
            6, Measurement.objects.filter(label="Lum_bc_Nref_median").count()
        )

    def test_unknown_measurement_and_method_are_refused(self):
        self.login()

        unknown_label = self.correct(
            {"label": "Fluo", "reference_type": "Nref", "method": "median"}
        )
        unknown_method = self.correct(
            {"label": "Lum", "reference_type": "Nref", "method": "max"}
        )

        self.assertEqual(400, unknown_label.status_code)
        self.assertIn('no measurement "Fluo"', unknown_label.json()[0])
        self.assertEqual(400, unknown_method.status_code)
        self.assertIn("method", unknown_method.json())

    def test_unknown_plate_is_not_found(self):
        self.login()
        url = reverse("correct_plate_background", args=[self.plate.id + 1000])

        response = self.client.post(
            url,
            {"label": "Lum", "reference_type": "Nref", "method": "median"},
            content_type="application/json",
        )

        self.assertEqual(404, response.status_code)
