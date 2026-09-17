"""
Tests for BaseMapper: going through the files, and creating plates and
measurement assignments.
"""

import os
import shutil
import tempfile
from os.path import join
from unittest import mock

from django.core.management.base import CommandError
from django.test import TestCase, override_settings

from core.models import (
    BarcodeSpecification,
    Experiment,
    ExperimentDetail,
    Plate,
    PlateDetail,
    PlateDimension,
    Project,
    WellDetail,
)
from importer.mappers import BaseMapper


class RecordingMapper(BaseMapper):
    """
    Remembers what run() hands to parse() and map().

    parse_calls example: [(("open file", "/tmp/a.csv", "ascii", "hello\\n"), {"xml_file": False})]
    map_calls example:   [(["row"], {"xml_file": False, "filename": "/tmp/a.csv"})]
    """

    def __init__(self, parse_result):
        self.parse_result = parse_result
        self.parse_calls = []
        self.map_calls = []

    def parse(self, file, **kwargs):
        if isinstance(file, str):
            received = file
        else:
            received = ("open file", file.name, file.encoding, file.read())
        self.parse_calls.append((received, dict(kwargs)))
        return self.parse_result

    def map(self, data, **kwargs):
        self.map_calls.append((data, dict(kwargs)))


@mock.patch.object(ExperimentDetail, "refresh")
@mock.patch.object(WellDetail, "refresh")
@mock.patch.object(PlateDetail, "refresh")
@mock.patch("importer.mappers.base.message")
class BaseMapperRunTest(TestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def write(self, name, content="hello\n"):
        path = join(self.folder, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as file:
            file.write(content)
        return path

    def test_a_text_file_is_opened_and_parsed(self, message, *refreshes):
        path = self.write("a.csv")
        mapper = RecordingMapper(parse_result=["row"])

        mapper.run(join(self.folder, "*.csv"), room_name="room_1")

        self.assertEqual(
            [
                (
                    ("open file", path, "ascii", "hello\n"),
                    {"room_name": "room_1", "xml_file": False},
                )
            ],
            mapper.parse_calls,
        )
        self.assertEqual(
            [(["row"], {"room_name": "room_1", "xml_file": False, "filename": path})],
            mapper.map_calls,
        )
        self.assertEqual(
            [
                mock.call(f"Processing file {path}...", "info", "room_1"),
                mock.call("Refreshing materialized views...", "info", "room_1"),
            ],
            message.call_args_list,
        )
        for refresh in refreshes:
            refresh.assert_called_once_with(concurrently=True)

    def test_an_xml_file_is_marked_as_xml(self, message, *refreshes):
        path = self.write("a.xml")
        mapper = RecordingMapper(parse_result=["row"])

        mapper.run(join(self.folder, "*.xml"))

        self.assertEqual({"xml_file": True}, mapper.parse_calls[0][1])
        self.assertEqual(
            [(["row"], {"xml_file": True, "filename": path})], mapper.map_calls
        )

    def test_a_parse_result_with_extra_information_is_split(self, message, *refreshes):
        # The M1000 mapper returns (entries, extra information)
        path = self.write("20240610-121212_demo_1.asc")
        mapper = RecordingMapper(parse_result=(["row"], {"barcode": "demo_1"}))

        mapper.run(join(self.folder, "*.asc"))

        self.assertEqual(
            [(["row"], {"xml_file": False, "barcode": "demo_1", "filename": path})],
            mapper.map_calls,
        )

    def test_txt_and_xlsx_files_are_parsed_by_name(self, message, *refreshes):
        path = self.write("241014_125455_241008MP-1_1.txt")
        mapper = RecordingMapper(parse_result={"results": []})

        mapper.run(join(self.folder, "*.txt"), room_name="room_1")

        self.assertEqual([(path, {"room_name": "room_1"})], mapper.parse_calls)
        self.assertEqual(
            [({"results": []}, {"room_name": "room_1", "filename": path})],
            mapper.map_calls,
        )

    def test_without_files_a_warning_is_shown(self, message, *refreshes):
        mapper = RecordingMapper(parse_result=[])
        pattern = join(self.folder, "*.csv")

        mapper.run(pattern, room_name="room_1")

        self.assertEqual([], mapper.map_calls)
        self.assertEqual(
            [
                mock.call(f"No files found that match {pattern}.", "warning", "room_1"),
                mock.call("Refreshing materialized views...", "info", "room_1"),
            ],
            message.call_args_list,
        )
        for refresh in refreshes:
            refresh.assert_called_once_with(concurrently=True)

    def test_files_in_sub_folders_are_found(self, message, *refreshes):
        top = self.write("a.csv")
        nested = self.write(join("sub", "b.csv"))

        found = BaseMapper.get_files(join(self.folder, "**", "*.csv"))

        self.assertEqual(sorted([top, nested]), sorted(found))


@mock.patch("importer.mappers.base.message")
class BaseMapperPlateTest(TestCase):
    fixtures = ["plate_dimensions"]

    def setUp(self):
        project = Project.objects.create(name="Project")
        self.experiment = Experiment.objects.create(name="Experiment", project=project)

    def test_the_experiment_comes_from_an_existing_barcode_specification(self, message):
        BarcodeSpecification.objects.create(
            prefix="2026Wagner12", experiment=self.experiment
        )

        plate = BaseMapper().create_plate_by_name_and_barcode(
            "Greiner_384PS_781904",
            "",
            "2026Wagner12",
            "384LDV_DMSO",
            experiment_name="Not used",
            room_name="room_1",
        )

        self.assertEqual("2026Wagner12", plate.barcode)
        self.assertEqual(self.experiment, plate.experiment)
        self.assertEqual("dim_384_16x24", plate.dimension.name)
        self.assertEqual(1, BarcodeSpecification.objects.count())
        message.assert_not_called()

    def test_a_missing_barcode_specification_is_created_for_the_experiment(
        self, message
    ):
        plate = BaseMapper().create_plate_by_name_and_barcode(
            "Corning_96",
            "",
            "ABC_1",
            "Source",
            experiment_name="Experiment",
            room_name="room_1",
        )

        specification = BarcodeSpecification.objects.get()
        self.assertEqual("ABC", specification.prefix)
        self.assertEqual(["North"], specification.sides)
        self.assertEqual(4, specification.number_of_plates)
        self.assertEqual(self.experiment, specification.experiment)
        self.assertEqual(self.experiment, plate.experiment)
        self.assertEqual("dim_96_8x12", plate.dimension.name)
        message.assert_called_once_with(
            "No barcode specification found for ABC_1. Creating it.",
            "warning",
            "room_1",
        )

    def test_without_barcode_specification_and_experiment_no_plate_is_created(
        self, message
    ):
        with self.assertRaises(CommandError) as raised:
            BaseMapper().create_plate_by_name_and_barcode(
                "Corning_96", "", "ABC_1", "Source", room_name="room_1"
            )

        self.assertEqual(
            "No barcode specification found for ABC_1 and no experiment name is provided. "
            "Please provide the experiment name in order to create the missing barcode "
            "specifications.",
            str(raised.exception),
        )
        # The map command shows the error, so it is not sent here too
        message.assert_not_called()
        self.assertFalse(Plate.objects.exists())

    def test_a_name_without_plate_size_is_refused(self, message):
        with self.assertRaises(CommandError) as raised:
            BaseMapper().get_plate_dimension("Plate_X", "", "Source_Y")

        self.assertEqual(
            "Could not determine plate dimensions for Plate_X: "
            "Cannot determine plate dimension from name:  Plate_X Source_Y.",
            str(raised.exception),
        )
        message.assert_not_called()

    def test_a_plate_size_without_plate_dimension_is_refused(self, message):
        PlateDimension.objects.filter(name="dim_96_8x12").delete()

        with self.assertRaises(CommandError) as raised:
            BaseMapper().get_plate_dimension("Corning_96", "", "Source")

        self.assertEqual(
            "No plate dimension found for Corning_96", str(raised.exception)
        )
        message.assert_not_called()


class BaseMapperMeasurementAssignmentTest(TestCase):
    fixtures = ["plate_dimensions"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def test_a_successful_assignment_with_the_file_is_created(self):
        path = join(self.folder, "20240610-121212_demo_1.asc")
        with open(path, "w") as file:
            file.write("A1\tSM1_1\t15\n")
        plate = Plate.objects.create(
            barcode="demo_1", dimension=PlateDimension.objects.get(name="dim_384_16x24")
        )

        with override_settings(MEDIA_ROOT=join(self.folder, "media")):
            assignment = BaseMapper().create_measurement_assignment(plate, path)

            self.assertEqual("success", assignment.status)
            self.assertEqual(plate, assignment.plate)
            self.assertEqual(path, assignment.filename)
            self.assertEqual(
                "20240610-121212_demo_1.asc", assignment.measurement_file.name
            )
            with assignment.measurement_file.open("r") as stored:
                self.assertEqual("A1\tSM1_1\t15\n", stored.read())
