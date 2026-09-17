"""
Tests for EchoMapper.map: turning parsed Echo transfers into plate mappings.

Plate.map itself is replaced by a recorder, so these tests check what the Echo
mapper hands over to it, not the mapping of wells.
"""

import shutil
import tempfile
from os.path import join
from unittest import mock

from django.test import TestCase, override_settings

from core.models import (
    Experiment,
    Plate,
    PlateDetail,
    PlateDimension,
    PlateMapping,
    Project,
    Well,
    WellDetail,
)
from importer.mappers import EchoMapper


def transfer(
    source_well, destination_well, source="SRC_A", destination="DST_1", **changes
):
    """
    One parsed Echo transfer, as EchoMapper.parse returns it.
    Example: transfer("A3", "B4", actual_volume="5") -> {"source_well": "A3", ...}
    """
    entry = {
        "source_plate_name": "384LDV_DMSO",
        "source_plate_barcode": source,
        "source_plate_type": "384LDV_DMSO",
        "source_well": source_well,
        "destination_plate_name": "Greiner_384PS_781904",
        "destination_plate_barcode": destination,
        "destination_well": destination_well,
        "actual_volume": "10",
        "transfer_status": "",
        "current_fluid_volume": "10.184",
        "DMSO": "98.744",
    }
    entry.update(changes)
    return entry


def describe(mapping):
    return (
        mapping.from_pos,
        mapping.to_pos,
        mapping.amount,
        mapping.status,
        mapping.map_type,
        mapping.current_amount,
        mapping.current_dmso,
    )


@mock.patch.object(WellDetail, "refresh")
@mock.patch.object(PlateDetail, "refresh")
@mock.patch("importer.mappers.echo.message")
class EchoMapTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.media = override_settings(MEDIA_ROOT=join(self.folder, "media"))
        self.media.enable()
        self.filename = join(self.folder, "Echo_Transfer-Echo_01_1.csv")
        with open(self.filename, "w") as file:
            file.write("transfer file\n")

        project = Project.objects.create(name="Project")
        Experiment.objects.create(name="Experiment", project=project)
        self.dimension = PlateDimension.objects.get(name="dim_384_16x24")
        for barcode in ("SRC_A", "SRC_B", "DST_1", "DST_2"):
            Plate.objects.create(barcode=barcode, dimension=self.dimension)
        # The source plates have all their wells, except where a test removes one
        for barcode in ("SRC_A", "SRC_B"):
            plate = Plate.objects.get(barcode=barcode)
            Well.objects.bulk_create(
                Well(plate=plate, position=position)
                for position in range(self.dimension.num_wells)
            )

        # Plate.map is replaced; every call is recorded as
        # (source barcode, target barcode, [described mappings])
        self.plate_map_calls = []
        patcher = mock.patch.object(
            Plate, "map", autospec=True, side_effect=self.record
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def tearDown(self):
        self.media.disable()
        shutil.rmtree(self.folder)

    def record(self, source_plate, mapping_list, target):
        self.plate_map_calls.append(
            (
                source_plate.barcode,
                target.barcode,
                [describe(item) for item in mapping_list],
            )
        )
        return True

    def run_map(self, data):
        EchoMapper().map(
            data,
            filename=self.filename,
            experiment_name="Experiment",
            room_name="room_1",
        )

    def position(self, well):
        return self.dimension.position(well)

    def test_the_transfers_of_one_plate_pair_are_mapped_together(
        self, message, plate_refresh, well_refresh
    ):
        self.run_map(
            [
                transfer("A3", "A3", DMSO="98.7%"),
                transfer(
                    "B4", "C5", actual_volume="2.5", current_fluid_volume="", DMSO=""
                ),
            ]
        )

        self.assertEqual(
            [
                (
                    "SRC_A",
                    "DST_1",
                    [
                        (
                            self.position("A3"),
                            self.position("A3"),
                            10.0,
                            "",
                            False,
                            10.184,
                            98.7,
                        ),
                        (
                            self.position("B4"),
                            self.position("C5"),
                            2.5,
                            "",
                            False,
                            None,
                            None,
                        ),
                    ],
                )
            ],
            self.plate_map_calls,
        )
        mapping = PlateMapping.objects.get()
        self.assertEqual(
            ("SRC_A", "DST_1"),
            (mapping.source_plate.barcode, mapping.target_plate.barcode),
        )
        self.assertEqual("Echo_Transfer-Echo_01_1.csv", mapping.mapping_file.name)
        self.assertEqual(
            [
                mock.call("Mapping SRC_A -> DST_1", "info", "room_1"),
                mock.call("Mapped SRC_A -> DST_1", "info", "room_1"),
            ],
            message.call_args_list,
        )
        # The views are refreshed once at the end of BaseMapper.run
        plate_refresh.assert_not_called()
        well_refresh.assert_not_called()

    def test_every_plate_pair_is_mapped_in_the_order_it_appears(
        self, message, *refreshes
    ):
        self.run_map(
            [
                transfer("A3", "A3"),
                transfer("A4", "A4", source="SRC_B"),
                transfer("A5", "A5", destination="DST_2"),
                transfer("A6", "A6"),
            ]
        )

        self.assertEqual(
            [("SRC_A", "DST_1"), ("SRC_B", "DST_1"), ("SRC_A", "DST_2")],
            [(source, target) for source, target, _ in self.plate_map_calls],
        )
        self.assertEqual(2, len(self.plate_map_calls[0][2]))
        self.assertEqual(3, PlateMapping.objects.count())

    def test_a_control_plate_as_source_maps_the_well_type(self, message, *refreshes):
        Plate.objects.filter(barcode="SRC_A").update(is_control_plate=True)

        self.run_map([transfer("A3", "A3")])

        self.assertTrue(self.plate_map_calls[0][2][0][4])

    def test_a_missing_destination_plate_is_created_once(self, message, *refreshes):
        self.run_map(
            [
                transfer("A3", "A3", destination="NEW_1"),
                transfer("A4", "A4", destination="NEW_1"),
            ]
        )

        plate = Plate.objects.get(barcode="NEW_1")
        self.assertEqual("dim_384_16x24", plate.dimension.name)
        self.assertEqual("Experiment", plate.experiment.name)
        self.assertEqual(
            [mock.call("Creating destination plate Greiner_384PS_781904, ")],
            [item for item in message.call_args_list if "Creating" in item.args[0]],
        )
        self.assertEqual(
            [("SRC_A", "NEW_1")], [(s, t) for s, t, _ in self.plate_map_calls]
        )

    def test_the_destination_plate_type_is_part_of_the_creation_message(
        self, message, *refreshes
    ):
        self.run_map(
            [
                transfer(
                    "A3", "A3", destination="NEW_1", destination_plate_type="Greiner"
                )
            ]
        )

        message.assert_any_call(
            "Creating destination plate Greiner_384PS_781904, Greiner"
        )

    def test_a_missing_source_plate_is_tried_again_and_then_reported_once(
        self, message, *refreshes
    ):
        self.run_map(
            [
                transfer("A3", "A3", source="MISSING"),
                transfer("A5", "A5", source="MISSING"),
                transfer("A4", "A4"),
            ]
        )

        self.assertEqual(
            [
                mock.call("Mapping SRC_A -> DST_1", "info", "room_1"),
                mock.call("Mapped SRC_A -> DST_1", "info", "room_1"),
                mock.call(
                    "2 transfers were not mapped, because these source plates "
                    "do not exist: MISSING",
                    "warning",
                    "room_1",
                ),
            ],
            message.call_args_list,
        )
        self.assertEqual(
            [("SRC_A", "DST_1")], [(s, t) for s, t, _ in self.plate_map_calls]
        )
        self.assertEqual(1, PlateMapping.objects.count())

    def test_the_same_report_is_not_mapped_again_onto_the_same_plates(
        self, message, *refreshes
    ):
        self.run_map([transfer("A3", "A3")])
        first_mapping = PlateMapping.objects.get()
        message.reset_mock()

        self.run_map([transfer("A3", "A3"), transfer("A3", "A3", destination="DST_2")])

        self.assertEqual(
            [
                mock.call("Mapping SRC_A -> DST_1", "info", "room_1"),
                mock.call(
                    "SRC_A -> DST_1: this report was already mapped on "
                    f"{first_mapping.created_at:%d.%m.%Y %H:%M}, so it was not mapped "
                    "again. To map it anew, delete the destination plate first.",
                    "error",
                    "room_1",
                ),
                # A plate pair that was not mapped before is still mapped
                mock.call("Mapping SRC_A -> DST_2", "info", "room_1"),
                mock.call("Mapped SRC_A -> DST_2", "info", "room_1"),
            ],
            message.call_args_list,
        )
        self.assertEqual(
            [("SRC_A", "DST_1"), ("SRC_A", "DST_2")],
            [(source, target) for source, target, _ in self.plate_map_calls],
        )
        self.assertEqual(2, PlateMapping.objects.count())

    def test_another_report_for_the_same_plates_is_not_reported(
        self, message, *refreshes
    ):
        self.run_map([transfer("A3", "A3")])
        with open(self.filename, "w") as file:
            file.write("another transfer file\n")
        message.reset_mock()

        self.run_map([transfer("A4", "A4")])

        self.assertEqual(
            [
                mock.call("Mapping SRC_A -> DST_1", "info", "room_1"),
                mock.call("Mapped SRC_A -> DST_1", "info", "room_1"),
            ],
            message.call_args_list,
        )

    def test_the_same_report_after_deleting_the_destination_plate_is_not_reported(
        self, message, *refreshes
    ):
        self.run_map([transfer("A3", "A3", destination="NEW_1")])
        Plate.objects.get(barcode="NEW_1").delete()
        message.reset_mock()

        self.run_map([transfer("A3", "A3", destination="NEW_1")])

        texts = [call.args[0] for call in message.call_args_list]
        self.assertFalse([text for text in texts if "already mapped" in text])
        self.assertEqual(1, PlateMapping.objects.count())

    def test_transfers_from_missing_source_wells_are_reported_once_per_plate_pair(
        self, message, *refreshes
    ):
        Well.objects.filter(
            plate__barcode="SRC_A",
            position__in=[self.position("L11"), self.position("A3")],
        ).delete()

        self.run_map(
            [
                transfer("L11", "A1"),
                transfer("A3", "A2"),
                transfer("L11", "A3"),
                transfer("A4", "A4"),
            ]
        )

        self.assertEqual(
            [
                mock.call("Mapping SRC_A -> DST_1", "info", "room_1"),
                mock.call("Mapped SRC_A -> DST_1", "info", "room_1"),
                mock.call(
                    "SRC_A -> DST_1: 3 transfers were not mapped, because these "
                    "source wells do not exist in SRC_A: A3, L11",
                    "warning",
                    "room_1",
                ),
            ],
            message.call_args_list,
        )
        # All transfers are still handed to Plate.map, which skips the missing wells
        self.assertEqual(4, len(self.plate_map_calls[0][2]))
        self.assertEqual(1, PlateMapping.objects.count())
