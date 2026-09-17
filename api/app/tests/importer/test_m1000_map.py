"""
Tests for M1000Mapper.map: writing parsed M1000 values as measurements.
"""

import shutil
import tempfile
from datetime import datetime
from os.path import join
from unittest import mock

from django.core.management.base import CommandError
from django.test import TestCase, override_settings

from core.models import (
    BarcodeSpecification,
    Experiment,
    Measurement,
    MeasurementAssignment,
    Plate,
    PlateDimension,
    Project,
    Well,
)
from importer.mappers import M1000Mapper

MEASURED_AT = datetime(2011, 11, 11, 11, 11, 11)


class M1000MapTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        media = override_settings(MEDIA_ROOT=join(self.folder, "media"))
        media.enable()
        self.addCleanup(media.disable)
        self.filename = join(self.folder, "20240610-121212_demo_1.asc")
        with open(self.filename, "w") as file:
            file.write("A1\tSM1_1\t15\n")

        project = Project.objects.create(name="Project")
        self.experiment = Experiment.objects.create(name="Experiment", project=project)
        self.dimension = PlateDimension.objects.get(name="dim_96_8x12")

        # One mock for the messages of m1000.py and base.py, so the order of all
        # messages is recorded, wherever the code sends them from.
        self.message = mock.Mock()
        for target in (
            "importer.mappers.m1000.message",
            "importer.mappers.base.message",
        ):
            patcher = mock.patch(target, self.message)
            patcher.start()
            self.addCleanup(patcher.stop)

    def tearDown(self):
        shutil.rmtree(self.folder)

    def run_map(self, entries, meta_data=None, **changes):
        data = {
            "barcode": "demo_1",
            "measurement_date": MEASURED_AT,
            "plate_description": None,
            "meta_data": meta_data or [{"Label": "Label1"}],
            "entries": entries,
        }
        kwargs = {
            "filename": self.filename,
            "experiment_name": "Experiment",
            "room_name": "room_1",
        }
        kwargs.update(changes)
        M1000Mapper().map(data, **kwargs)

    def measurements(self):
        """[(well position, label, value, identifier), ...] sorted."""
        return sorted(
            Measurement.objects.values_list(
                "well__position", "label", "value", "identifier"
            )
        )

    def test_the_values_become_measurements_of_the_wells(self):
        plate = Plate.objects.create(barcode="demo_1", dimension=self.dimension)
        existing_well = Well.objects.create(plate=plate, position=0)

        self.run_map(
            [
                {"position": "A1", "identifier": "SM1_1", "values": [15.0]},
                {"position": "B2", "identifier": "SM1_14", "values": [999.0]},
            ]
        )

        self.assertEqual(
            [(0, "Label1", 15.0, "SM1_1"), (13, "Label1", 999.0, "SM1_14")],
            self.measurements(),
        )
        self.assertEqual(2, plate.wells.count())
        self.assertEqual(existing_well, Measurement.objects.get(well__position=0).well)

        assignment = MeasurementAssignment.objects.get()
        self.assertEqual(
            (plate, "success", self.filename),
            (assignment.plate, assignment.status, assignment.filename),
        )
        for measurement in Measurement.objects.all():
            self.assertEqual(assignment, measurement.measurement_assignment)
            self.assertEqual(MEASURED_AT, measurement.measured_at)
        self.message.assert_not_called()

    def test_the_labels_come_from_the_meta_data_in_order(self):
        Plate.objects.create(barcode="demo_1", dimension=self.dimension)

        self.run_map(
            [{"position": "A1", "identifier": "SM1_1", "values": [1.0, 2.0]}],
            meta_data=[{"Label": "First"}, {"Label": "Second"}],
        )

        self.assertEqual(
            [(0, "First", 1.0, "SM1_1"), (0, "Second", 2.0, "SM1_1")],
            self.measurements(),
        )

    def test_given_measurement_names_replace_the_labels(self):
        Plate.objects.create(barcode="demo_1", dimension=self.dimension)

        self.run_map(
            [{"position": "A1", "identifier": "SM1_1", "values": [1.0, 2.0]}],
            measurement_name="Lum,Fluo",
        )

        self.assertEqual(
            [(0, "Fluo", 2.0, "SM1_1"), (0, "Lum", 1.0, "SM1_1")],
            self.measurements(),
        )

    def test_fewer_measurement_names_than_value_columns_are_refused(self):
        Plate.objects.create(barcode="demo_1", dimension=self.dimension)

        with self.assertRaisesMessage(
            CommandError, "2 value columns, but only 1 measurement names"
        ):
            self.run_map(
                [{"position": "A1", "identifier": "SM1_1", "values": [1.0, 2.0]}],
                measurement_name="Lum",
            )

        self.assertFalse(Measurement.objects.exists())

    def test_a_footer_without_a_label_is_refused(self):
        Plate.objects.create(barcode="demo_1", dimension=self.dimension)

        with self.assertRaisesMessage(
            CommandError, "its footer does not name every label"
        ):
            self.run_map(
                [{"position": "A1", "identifier": "SM1_1", "values": [1.0]}],
                meta_data=[{"Integration time": "1000 ms"}],
            )

    def test_a_missing_plate_is_created_for_the_experiment(self):
        self.run_map(
            [
                {"position": "A1", "identifier": "SM1_1", "values": [15.0]},
                {"position": "A2", "identifier": "SM1_2", "values": [16.0]},
            ]
        )

        plate = Plate.objects.get(barcode="demo_1")
        self.assertEqual("dim_96_8x12", plate.dimension.name)
        self.assertEqual(self.experiment, plate.experiment)
        specification = BarcodeSpecification.objects.get()
        self.assertEqual(
            ("demo", ["North"], 4, self.experiment),
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
                    "Plate with barcode demo_1 does not exist. Creating it.",
                    "warning",
                    "room_1",
                )
            ],
            self.message.call_args_list,
        )
        self.assertEqual(2, Measurement.objects.count())

    def test_mapping_the_same_values_again_keeps_one_measurement_per_label(self):
        Plate.objects.create(barcode="demo_1", dimension=self.dimension)
        data = [{"position": "A1", "identifier": "SM1_1", "values": [15.0]}]

        self.run_map(data)
        self.run_map(data)

        self.assertEqual([(0, "Label1", 15.0, "SM1_1")], self.measurements())
