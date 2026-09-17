"""
Tests for the command that adds supplier data from a CSV file to the compounds.
"""

import shutil
import tempfile
from os.path import join

from django.core.management import call_command
from django.test import TestCase

from compoundlib.models import Compound


class UpdateCompoundsTest(TestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)

    def test_a_file_saved_by_excel_is_read(self):
        # Excel writes an invisible BOM character before the first column name
        compound = Compound.objects.create(name="Aspirin", data={"Plate": "P1"})
        path = join(self.folder, "compounds.tsv")
        with open(path, "w", encoding="utf-8") as file:
            file.write("\ufeffCompoundName\tTarget\nAspirin\tCOX\n")

        call_command("update_compounds", input_file=path)

        compound.refresh_from_db()
        self.assertEqual(
            {"Plate": "P1", "CompoundName": "Aspirin", "Target": "COX"}, compound.data
        )
