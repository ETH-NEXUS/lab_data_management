"""
Tests for MicroscopeMapper.map: writing parsed C10 results as measurements.
"""

import shutil
import tempfile
from datetime import datetime
from os.path import join
from unittest import mock

from django.test import TestCase, override_settings

from core.models import (
    BarcodeSpecification,
    Experiment,
    Measurement,
    MeasurementAssignment,
    Plate,
    PlateDimension,
    Project,
)
from importer.mappers import MicroscopeMapper


class MicroscopeMapTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        media = override_settings(MEDIA_ROOT=join(self.folder, "media"))
        media.enable()
        self.addCleanup(media.disable)
        self.filename = join(self.folder, "241014_125455_241008MP-1_1.txt")
        with open(self.filename, "w") as file:
            file.write("Results\n")

        project = Project.objects.create(name="Project")
        self.experiment = Experiment.objects.create(name="Experiment", project=project)
        self.dimension = PlateDimension.objects.get(name="dim_384_16x24")

        # One mock for the messages of microscope.py and base.py, so all
        # messages are recorded in order, wherever the code sends them from.
        self.message = mock.Mock()
        for target in (
            "importer.mappers.microscope.message",
            "importer.mappers.base.message",
        ):
            patcher = mock.patch(target, self.message)
            patcher.start()
            self.addCleanup(patcher.stop)

    def tearDown(self):
        shutil.rmtree(self.folder)

    def run_map(
        self, results, layout=None, barcode="241008MP-1_1", date="241014", time="125455"
    ):
        data = {
            "metadata": {},
            "results": results,
            "date": date,
            "time": time,
            "barcode": barcode,
            "layout": layout if layout is not None else {},
        }
        MicroscopeMapper().map(
            data,
            filename=self.filename,
            experiment_name="Experiment",
            room_name="room_1",
        )

    def measurements(self):
        """[(well position, label, value), ...] sorted."""
        return sorted(
            Measurement.objects.values_list("well__position", "label", "value")
        )

    def test_the_numbers_become_measurements_of_the_wells(self):
        plate = Plate.objects.create(barcode="241008MP-1_1", dimension=self.dimension)

        self.run_map(
            [
                {"Well": "A1", "Lum": "16727"},
                {"Well": "A2", "Lum": "1.5E+03"},
                {"Well": "A3", "Lum": "n/a"},
                {"Well ID": "SPL1", "Well": "B1", "Lum": 16727, "Note": "text"},
            ]
        )

        self.assertEqual(
            [(0, "Lum", 16727.0), (1, "Lum", 1500.0), (24, "Lum", 16727.0)],
            self.measurements(),
        )
        # A well is created even when none of its values is a number
        self.assertEqual(
            [0, 1, 2, 24], sorted(plate.wells.values_list("position", flat=True))
        )
        for measurement in Measurement.objects.all():
            self.assertEqual(
                datetime(2024, 10, 14, 12, 54, 55), measurement.measured_at
            )
        assignment = MeasurementAssignment.objects.get()
        self.assertEqual(
            (plate, "success", self.filename),
            (assignment.plate, assignment.status, assignment.filename),
        )
        self.message.assert_not_called()

    def test_empty_rows_and_header_rows_are_not_wells(self):
        plate = Plate.objects.create(barcode="241008MP-1_1", dimension=self.dimension)

        self.run_map(
            [
                {"Well": None, "Lum": "1"},
                {"Well": "", "Lum": "2"},
                {"Well": "Well", "Lum": "Lum"},
            ]
        )

        self.assertFalse(plate.wells.exists())
        self.assertFalse(Measurement.objects.exists())
        self.assertEqual(1, MeasurementAssignment.objects.count())

    def test_the_layout_sets_the_well_types(self):
        plate = Plate.objects.create(barcode="241008MP-1_1", dimension=self.dimension)

        self.run_map(
            [
                {"Well": "A1", "Lum": "1"},
                {"Well": "B1", "Lum": "2"},
                {"Well": "A2", "Lum": "3"},
            ],
            layout={"A1": "P", "B1": "N"},
        )

        well_types = dict(plate.wells.values_list("position", "type__name"))
        self.assertEqual({0: "P", 24: "N", 1: "C"}, well_types)

    def test_a_missing_plate_is_created_for_the_experiment(self):
        self.run_map(
            [
                {"Well": "A1", "Lum": "1"},
                {"Well": "A2", "Lum": "2"},
                {"Well": "A3", "Lum": "3"},
            ]
        )

        plate = Plate.objects.get(barcode="241008MP-1_1")
        self.assertEqual("dim_96_8x12", plate.dimension.name)
        self.assertEqual(self.experiment, plate.experiment)
        specification = BarcodeSpecification.objects.get()
        self.assertEqual(
            ("241008MP-1", ["North"], 4, self.experiment),
            (
                specification.prefix,
                specification.sides,
                specification.number_of_plates,
                specification.experiment,
            ),
        )
        self.assertEqual(
            [
                mock.call(
                    "Plate with barcode 241008MP-1_1 does not exist. Creating it.",
                    "warning",
                    "room_1",
                )
            ],
            self.message.call_args_list,
        )
        self.assertEqual(3, Measurement.objects.count())

    def test_an_unknown_date_stops_before_anything_is_stored(self):
        with self.assertRaisesMessage(
            ValueError,
            f"Cannot read the measurement date of {self.filename}: "
            "date '14.10.2024', time '12:45:28'.",
        ):
            self.run_map(
                [{"Well": "A1", "Lum": "1"}], date="14.10.2024", time="12:45:28"
            )

        self.assertFalse(Plate.objects.exists())
        self.assertFalse(Measurement.objects.exists())
        self.assertFalse(MeasurementAssignment.objects.exists())

    def test_mapping_the_same_results_again_keeps_one_measurement_per_label(self):
        Plate.objects.create(barcode="241008MP-1_1", dimension=self.dimension)
        results = [{"Well": "A1", "Lum": "16727"}]

        self.run_map(results)
        self.run_map(results)

        self.assertEqual([(0, "Lum", 16727.0)], self.measurements())
