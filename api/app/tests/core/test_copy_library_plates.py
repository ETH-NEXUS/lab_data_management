"""
Tests for copying library plates, as the "Copy selected plates" admin action does.
"""

from unittest.mock import patch

from django.test import TestCase

from compoundlib.models import CompoundLibrary
from core.models import Plate, PlateDetail, PlateDimension, Well, WellDetail
from core.utils import copy_library_plates


@patch.object(WellDetail, "refresh")
@patch.object(PlateDetail, "refresh")
class CopyLibraryPlatesTest(TestCase):
    fixtures = ("well_types",)

    def setUp(self):
        dimension = PlateDimension.objects.create(name="dim_3x2", cols=3, rows=2)
        library = CompoundLibrary.objects.create(name="Test library")
        self.source = Plate.objects.create(
            barcode="OLD_PLATE", dimension=dimension, library=library, archived=True
        )
        Well.objects.create(plate=self.source, position=0)

    def test_the_copy_of_an_archived_plate_is_not_archived(
        self, plate_refresh, well_refresh
    ):
        copy = copy_library_plates([self.source], 500)[0]
        copy.refresh_from_db()
        self.assertFalse(copy.archived)

    def test_the_source_plate_stays_archived(self, plate_refresh, well_refresh):
        copy_library_plates([self.source], 500)
        self.source.refresh_from_db()
        self.assertTrue(self.source.archived)
