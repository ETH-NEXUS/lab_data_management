"""
Tests for marking source wells that are running low while a plate is mapped.

The Echo reports the fill level of a source well with every transfer, also
when the transfer itself failed. Mapping has to act on that report right away,
including the very first transfer out of a well.
"""

from django.test import TestCase

from compoundlib.models import Compound, CompoundLibrary
from core.mapping import Mapping, MappingList
from core.models import Plate, PlateDimension, Threshold, Well, WellCompound


class PlateMapThresholdTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        Threshold.objects.create(amount=2.5, dmso=80)
        self.dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        self.library = CompoundLibrary.objects.create(name="Test library")
        self.source_plate = Plate.objects.create(
            barcode="SOURCE", dimension=self.dimension, library=self.library
        )
        self.target_plate = Plate.objects.create(
            barcode="TARGET", dimension=self.dimension
        )
        self.compound = Compound.objects.create(name="Test compound")

    def map_well(self, position, current_amount, current_dmso, amount=30):
        """
        Maps one well and returns it, as the Echo import would after reading a
        row like {"current_fluid_volume": 2.0, "DMSO": 95}.
        """
        well = Well.objects.create(plate=self.source_plate, position=position)
        WellCompound.objects.create(well=well, compound=self.compound, amount=1000)

        mapping_list = MappingList(target=self.target_plate)
        mapping_list.add(
            Mapping(
                from_pos=position,
                to_pos=position,
                amount=amount,
                current_amount=current_amount,
                current_dmso=current_dmso,
            )
        )
        self.source_plate.map(mapping_list, self.target_plate)

        well.refresh_from_db()
        self.source_plate.refresh_from_db()
        return well

    def test_a_low_volume_is_marked_on_the_first_transfer(self):
        well = self.map_well(0, current_amount=2.0, current_dmso=95)
        self.assertEqual("empty", well.status)
        self.assertEqual("empty_wells", self.source_plate.status)

    def test_a_failed_transfer_reporting_zeros_is_marked(self):
        well = self.map_well(1, current_amount=0, current_dmso=0)
        self.assertEqual("empty", well.status)

    def test_a_low_dmso_is_marked(self):
        well = self.map_well(2, current_amount=9.5, current_dmso=70)
        self.assertEqual("empty", well.status)

    def test_a_well_that_is_fine_is_not_marked(self):
        well = self.map_well(3, current_amount=9.5, current_dmso=95)
        self.assertIsNone(well.status)
        self.assertIsNone(self.source_plate.status)

    def test_a_mapping_without_any_reading_marks_nothing(self):
        """A plate copy or a csv mapping carries no fill level at all."""
        well = Well.objects.create(plate=self.source_plate, position=4)
        WellCompound.objects.create(well=well, compound=self.compound, amount=1000)

        self.source_plate.copy(self.target_plate, 30)

        well.refresh_from_db()
        self.source_plate.refresh_from_db()
        self.assertIsNone(well.status)
        self.assertIsNone(self.source_plate.status)

    def test_a_low_volume_is_marked_on_a_repeated_transfer_too(self):
        well = self.map_well(5, current_amount=9.5, current_dmso=95)
        self.assertIsNone(well.status)

        mapping_list = MappingList(target=self.target_plate)
        mapping_list.add(
            Mapping(
                from_pos=5, to_pos=5, amount=30, current_amount=2.0, current_dmso=95
            )
        )
        self.source_plate.map(mapping_list, self.target_plate)

        well.refresh_from_db()
        self.assertEqual("empty", well.status)

    def test_wells_of_a_plate_without_a_library_are_left_alone(self):
        plate = Plate.objects.create(barcode="NO_LIBRARY", dimension=self.dimension)
        well = Well.objects.create(plate=plate, position=0)
        WellCompound.objects.create(well=well, compound=self.compound, amount=1000)

        mapping_list = MappingList(target=self.target_plate)
        mapping_list.add(
            Mapping(from_pos=0, to_pos=0, amount=30, current_amount=0, current_dmso=0)
        )
        plate.map(mapping_list, self.target_plate)

        well.refresh_from_db()
        plate.refresh_from_db()
        self.assertIsNone(well.status)
        self.assertIsNone(plate.status)

    def test_a_status_set_by_hand_is_kept(self):
        # Set on the plate object: the Echo import also loads the plate with its status.
        self.source_plate.status = "disposed"
        self.source_plate.save()
        well = self.map_well(0, current_amount=0, current_dmso=0)
        self.assertEqual("empty", well.status)
        self.assertEqual("disposed", self.source_plate.status)
