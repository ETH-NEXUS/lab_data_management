"""
Tests for the command that removes the BOM character from compound names.
"""

from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from compoundlib.models import Compound
from core.models import Plate, PlateDimension, Well, WellCompound


class RemoveBomFromCompoundNamesTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        dimension = PlateDimension.objects.get(name="dim_96_8x12")
        plate = Plate.objects.create(barcode="CTRL_1", dimension=dimension)
        self.well_a1 = Well.objects.create(plate=plate, position=0)
        self.well_a2 = Well.objects.create(plate=plate, position=1)

    def run_command(self, *args):
        output = StringIO()
        call_command("remove_bom_from_compound_names", *args, stdout=output)
        return output.getvalue()

    def test_a_bom_compound_is_merged_into_the_compound_with_the_clean_name(self):
        dmso = Compound.objects.create(name="DMSO")
        bom_dmso = Compound.objects.create(name="\ufeffDMSO")
        WellCompound.objects.create(well=self.well_a1, compound=bom_dmso, amount=5)

        output = self.run_command()

        self.assertIn(f"Merged: '<BOM>DMSO' (id {bom_dmso.id}, 1 wells)", output)
        self.assertEqual(
            ["DMSO"], list(Compound.objects.values_list("name", flat=True))
        )
        self.assertEqual(
            [(self.well_a1.id, dmso.id, 5.0)],
            list(WellCompound.objects.values_list("well", "compound", "amount")),
        )

    def test_a_well_with_both_compounds_keeps_one_entry_with_both_amounts(self):
        dmso = Compound.objects.create(name="DMSO")
        bom_dmso = Compound.objects.create(name="\ufeffDMSO")
        WellCompound.objects.create(well=self.well_a1, compound=dmso, amount=2)
        WellCompound.objects.create(well=self.well_a1, compound=bom_dmso, amount=5)

        self.run_command()

        self.assertEqual(
            [(self.well_a1.id, dmso.id, 7.0)],
            list(WellCompound.objects.values_list("well", "compound", "amount")),
        )

    def test_a_bom_compound_without_a_clean_twin_is_renamed(self):
        bom_compound = Compound.objects.create(name="\ufeffCompound1")
        WellCompound.objects.create(well=self.well_a2, compound=bom_compound)

        output = self.run_command()

        self.assertIn("Renamed: '<BOM>Compound1'", output)
        bom_compound.refresh_from_db()
        self.assertEqual("Compound1", bom_compound.name)

    def test_a_dry_run_saves_nothing(self):
        Compound.objects.create(name="DMSO")
        bom_dmso = Compound.objects.create(name="\ufeffDMSO")
        WellCompound.objects.create(well=self.well_a1, compound=bom_dmso)

        output = self.run_command("--dry-run")

        self.assertIn("Merged: '<BOM>DMSO'", output)
        self.assertIn("Dry run: nothing was saved.", output)
        self.assertEqual(2, Compound.objects.count())
        self.assertEqual(bom_dmso.id, WellCompound.objects.get().compound_id)

    def test_without_bom_names_nothing_happens(self):
        Compound.objects.create(name="DMSO")

        output = self.run_command()

        self.assertEqual("No compound names with a BOM character.\n", output)
