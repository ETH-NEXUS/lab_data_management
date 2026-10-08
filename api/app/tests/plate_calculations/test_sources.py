"""
Where the measurements of a plate come from: the name of the file an imported
one was read from, nothing for a calculated one. Only for logged in users.
"""

from datetime import datetime

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from core.models import (
    Experiment,
    Measurement,
    MeasurementAssignment,
    Plate,
    PlateDimension,
    Project,
    Well,
    WellType,
)

READ = datetime(2026, 9, 30, 15, 46, 54)
FILE_PATH = "/data/projects/Run/M1000 Output/Evaluated/093026-154654_RKS_1.asc"


class SourcesViewTest(TestCase):
    def setUp(self):
        project = Project.objects.create(name="P1")
        experiment = Experiment.objects.create(name="Screen 1", project=project)
        dimension = PlateDimension.objects.create(name="dim_1x2", rows=1, cols=2)
        self.plate = Plate.objects.create(
            barcode="RKS_1", dimension=dimension, experiment=experiment
        )
        assignment = MeasurementAssignment.objects.create(
            plate=self.plate, filename=FILE_PATH
        )
        well_type = WellType.objects.create(name="N", description="negative")
        for position in range(2):
            well = Well.objects.create(
                plate=self.plate, position=position, type=well_type
            )
            # Imported, also the one the software of the instrument calculated
            for label in ["Lum1", " Lum1_corrected"]:
                Measurement.objects.create(
                    well=well,
                    label=label,
                    value=100.0,
                    measured_at=READ,
                    measurement_assignment=assignment,
                )
            Measurement.objects.create(
                well=well, label="Lum1_log10", value=2.0, measured_at=READ
            )

    def sources(self):
        url = reverse("plate_measurement_sources", args=[self.plate.id])
        return self.client.get(url)

    def test_without_login_nothing_is_shown(self):
        self.assertEqual(403, self.sources().status_code)

    def test_the_file_name_of_imported_measurements_and_none_for_the_others(self):
        self.client.force_login(User.objects.create_user("tester"))

        response = self.sources()

        self.assertEqual(200, response.status_code)
        self.assertEqual(
            {
                "Lum1": "093026-154654_RKS_1.asc",
                " Lum1_corrected": "093026-154654_RKS_1.asc",
                "Lum1_log10": None,
            },
            response.json(),
        )
