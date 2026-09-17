"""
Tests for the command that exports the objects of a model to a YAML file.
"""

import shutil
import tempfile
from os.path import join

import yaml
from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase

from compoundlib.models import Compound


class ExportCommandTest(TestCase):
    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.output_file = join(self.folder, "compounds.yaml")
        Compound.objects.create(name="Aspirin")
        Compound.objects.create(name="Caffeine")

    def tearDown(self):
        shutil.rmtree(self.folder)

    def exported_names(self):
        with open(self.output_file) as file:
            return [item["fields"]["name"] for item in yaml.safe_load(file)]

    def test_all_objects_are_exported(self):
        call_command("export", "compoundlib", "Compound", output_file=self.output_file)

        self.assertEqual(["Aspirin", "Caffeine"], sorted(self.exported_names()))

    def test_the_filters_select_the_objects(self):
        call_command(
            "export",
            "compoundlib",
            "Compound",
            "--filter",
            "name__startswith=Asp",
            "--filter",
            "name__endswith=rin",
            output_file=self.output_file,
        )

        self.assertEqual(["Aspirin"], self.exported_names())

    def test_a_filter_without_equals_sign_is_refused(self):
        with self.assertRaisesMessage(CommandError, "is not field=value"):
            call_command(
                "export",
                "compoundlib",
                "Compound",
                "--filter",
                "name",
                output_file=self.output_file,
            )
