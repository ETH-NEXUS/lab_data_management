"""
The plate information of "add experiment data": one row per plate and measurement
label. A plate measured with two labels keeps both rows when it is saved.
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
    PlateInfo,
    Project,
    Well,
    WellType,
)

MEASURED_AT = datetime(2025, 5, 16, 10, 0)


class PlateInfoTest(TestCase):
    def setUp(self):
        self.client.force_login(User.objects.create_user("tester"))
        project = Project.objects.create(name="P1")
        self.experiment = Experiment.objects.create(name="Screen 1", project=project)
        # The prefill looks at a well in the middle of the plate, so 2x2 at least
        dimension = PlateDimension.objects.create(name="dim_2x2", rows=2, cols=2)
        plate = Plate.objects.create(
            barcode="SP_1", dimension=dimension, experiment=self.experiment
        )
        well_type = WellType.objects.create(name="C", description="compound")
        for position in range(4):
            well = Well.objects.create(plate=plate, position=position, type=well_type)
            for label in ["Lum", "Fluo"]:
                Measurement.objects.create(
                    well=well, label=label, value=1, measured_at=MEASURED_AT
                )
        PlateDetail.refresh()

    def row(self, label, condition):
        return {
            "plate_barcode": "SP_1",
            "lib_plate_barcode": "LIB_1",
            "measurement_label": label,
            "measurement_timestamp": MEASURED_AT.isoformat(),
            "replicate": "1",
            "cell_type": "SW620",
            "condition": condition,
        }

    def save(self, rows):
        return self.client.post(
            reverse("save_plate_info"),
            {"experiment_id": self.experiment.id, "plate_info": rows},
            content_type="application/json",
        )

    def prefill(self):
        response = self.client.get(
            reverse("prefillPlateInfo"), {"experiment_id": self.experiment.id}
        )
        return response.json()["plate_info"]

    def test_both_labels_of_a_plate_are_saved(self):
        self.save([self.row("Lum", "KO"), self.row("Fluo", "WT")])

        saved = PlateInfo.objects.order_by("label").values_list("label", "condition")
        self.assertEqual([("Fluo", "WT"), ("Lum", "KO")], list(saved))

    def test_saving_again_changes_the_row_of_the_same_label(self):
        self.save([self.row("Lum", "KO")])
        self.save([self.row("Lum", "WT")])

        saved = PlateInfo.objects.values_list("label", "condition")
        self.assertEqual([("Lum", "WT")], list(saved))

    def test_the_form_shows_a_label_that_was_not_saved_yet(self):
        self.save([self.row("Lum", "KO")])

        rows = self.prefill()

        conditions = {row["measurement_label"]: row["condition"] for row in rows}
        self.assertEqual({"Lum": "KO", "Fluo": ""}, conditions)
