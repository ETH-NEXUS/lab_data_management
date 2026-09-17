"""
Tests that run a mapper over instrument files, from the files to the database,
the same way the `map` command does it (file patterns from ldm.yaml, real
Plate.map, real measurements).
"""

import shutil
import tempfile
from datetime import datetime
from os import makedirs
from os.path import dirname, join
from unittest import mock

from django.test import TestCase, override_settings

from compoundlib.models import Compound
from core.models import (
    BarcodeSpecification,
    Experiment,
    Measurement,
    MeasurementAssignment,
    Plate,
    PlateDimension,
    PlateMapping,
    Project,
    WellCompound,
    WellWithdrawal,
)
from importer.config import Config
from importer.mappers import EchoMapper, M1000Mapper

# Real M1000 files, as the reader writes them
INSTRUMENT_FILES = join(dirname(__file__), "instrument_files")
M1000_ONE_LABEL_FILE = join(INSTRUMENT_FILES, "20210902-131750_BAF210901_1.asc")
M1000_TWO_LABELS_FILE = join(INSTRUMENT_FILES, "20191205-101721_221212AK_1.asc")

# An Echo report from source plate LLD_4541_C (test fixture) to one destination plate
ECHO_REPORT = """Run ID,14618
Run Date/Time,02/09/2021 10:30:32
Application Name,
Application Version,1.0.0.0
Protocol Name,
User Name,cellario


[DETAILS]
Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Source Concentration,Source Concentration Units,Destination Plate Name,Destination Plate Barcode,Destination Well,Destination Concentration,Destination Concentration Units,Compound Name,Transfer Volume,Actual Volume,Transfer Status,Current Fluid Height,Current Fluid Volume,% DMSO
384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A10,0,N/A,Corning_384_3577,{barcode},A10,0,N/A,N/A,30,30,,2.003,9.959,99.203
384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A11,0,N/A,Corning_384_3577,{barcode},A11,0,N/A,N/A,30,30,,1.968,9.786,99.008
384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A12,0,N/A,Corning_384_3577,{barcode},A12,0,N/A,N/A,30,30,,2.001,9.959,99.348
"""


class MapperRunTest(TestCase):
    fixtures = ["plate_dimensions", "well_types", "test/compound_library"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        media = override_settings(MEDIA_ROOT=join(self.folder, "media"))
        media.enable()
        self.addCleanup(media.disable)
        self.project = Project.objects.create(name="Project")
        self.dimension = PlateDimension.objects.get(name="dim_384_16x24")

    def tearDown(self):
        shutil.rmtree(self.folder)

    def write(self, path, text):
        makedirs(dirname(path), exist_ok=True)
        with open(path, "w") as file:
            file.write(text)

    def test_echo_reports_in_sub_folders_are_mapped(self):
        experiment = Experiment.objects.create(name="Echo", project=self.project)
        source_plate = Plate.objects.get(barcode="LLD_4541_C")
        compound = Compound.objects.create(name="Compound A")
        source_positions = [
            self.dimension.position(well) for well in ("A10", "A11", "A12")
        ]
        for position in source_positions:
            well = source_plate.well_at(position, create_if_not_exist=True)
            WellCompound.objects.create(well=well, compound=compound)

        echo_folder = join(self.folder, "echo")
        for barcode in ("P1", "P2", "P3"):
            BarcodeSpecification.objects.create(prefix=barcode, experiment=experiment)
            self.write(
                join(echo_folder, barcode, "ID-123-transfer-Echo_01_123.csv"),
                ECHO_REPORT.format(barcode=barcode),
            )
            # Files that do not fit the file pattern are ignored
            self.write(join(echo_folder, barcode, "something.csv"), "")
            self.write(join(echo_folder, barcode, "dfg-transfer-file.tsv"), "")

        file_pattern = Config.current.importer.echo.default.file_blob
        EchoMapper().run(join(echo_folder, file_pattern))

        self.assertEqual(
            [("LLD_4541_C", "P1"), ("LLD_4541_C", "P2"), ("LLD_4541_C", "P3")],
            sorted(
                PlateMapping.objects.values_list(
                    "source_plate__barcode", "target_plate__barcode"
                )
            ),
        )
        for barcode in ("P1", "P2", "P3"):
            plate = Plate.objects.get(barcode=barcode)
            self.assertEqual(experiment, plate.experiment)
            self.assertEqual(
                [(position, "Compound A") for position in source_positions],
                sorted(
                    WellCompound.objects.filter(well__plate=plate).values_list(
                        "well__position", "compound__name"
                    )
                ),
            )
        # Every source well gave 30 nL to each of the three destination plates
        a10 = source_plate.well_at(source_positions[0])
        self.assertEqual(
            [(30.0, 9.959, 99.203)] * 3,
            list(
                WellWithdrawal.objects.filter(well=a10).values_list(
                    "amount", "current_amount", "current_dmso"
                )
            ),
        )

    def test_a_transfer_from_a_missing_source_well_is_skipped_and_reported(self):
        experiment = Experiment.objects.create(name="Echo", project=self.project)
        BarcodeSpecification.objects.create(prefix="P1", experiment=experiment)
        source_plate = Plate.objects.get(barcode="LLD_4541_C")
        compound = Compound.objects.create(name="Compound A")
        # A12 is not in the database
        for well_name in ("A10", "A11"):
            well = source_plate.well_at(
                self.dimension.position(well_name), create_if_not_exist=True
            )
            WellCompound.objects.create(well=well, compound=compound)
        source_plate.wells.filter(position=self.dimension.position("A12")).delete()
        path = join(self.folder, "ID-123-transfer-Echo_01_123.csv")
        self.write(path, ECHO_REPORT.format(barcode="P1"))

        with mock.patch("importer.mappers.echo.message") as message:
            EchoMapper().run(path, room_name="room_1")

        message.assert_any_call(
            "LLD_4541_C -> P1: 1 transfers were not mapped, because these source "
            "wells do not exist in LLD_4541_C: A12",
            "warning",
            "room_1",
        )
        plate = Plate.objects.get(barcode="P1")
        self.assertEqual(
            [self.dimension.position("A10"), self.dimension.position("A11")],
            sorted(
                WellCompound.objects.filter(well__plate=plate).values_list(
                    "well__position", flat=True
                )
            ),
        )

    def test_m1000_files_become_measurements(self):
        experiment = Experiment.objects.create(name="M1000", project=self.project)
        for number in (1, 2, 3):
            Plate.objects.create(
                barcode=f"BAF210901_{number}",
                dimension=self.dimension,
                experiment=experiment,
            )
            shutil.copy(
                M1000_ONE_LABEL_FILE,
                join(self.folder, f"20210902-131750_BAF210901_{number}.asc"),
            )

        file_pattern = Config.current.importer.m1000.default.file_blob
        M1000Mapper().run(join(self.folder, file_pattern))

        # 3 plates with 384 wells and one value per well
        self.assertEqual(1152, Measurement.objects.count())
        self.assertEqual(
            ["Label1"],
            list(Measurement.objects.values_list("label", flat=True).distinct()),
        )
        self.assertEqual(3, MeasurementAssignment.objects.count())
        a1 = Measurement.objects.get(
            well__plate__barcode="BAF210901_1", well__position=0
        )
        self.assertEqual(
            (11115.0, "NC1", datetime(2021, 9, 2, 13, 20, 1)),
            (a1.value, a1.identifier, a1.measured_at),
        )

    def test_the_labels_of_an_m1000_file_with_two_values_per_well(self):
        # The file has the columns "Acceptor" and "Donor", but its footer lists
        # the label "Donor" first: the labels are read in reverse order.
        plate = Plate.objects.create(barcode="221212AK_1", dimension=self.dimension)
        shutil.copy(M1000_TWO_LABELS_FILE, self.folder)

        file_pattern = Config.current.importer.m1000.default.file_blob
        M1000Mapper().run(join(self.folder, file_pattern))

        a1 = Measurement.objects.filter(well__plate=plate, well__position=0)
        self.assertEqual(
            {"Acceptor": 24672.0, "Donor": 7395.0},
            dict(a1.values_list("label", "value")),
        )
        self.assertEqual(2 * plate.wells.count(), Measurement.objects.count())
