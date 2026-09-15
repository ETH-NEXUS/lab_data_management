"""
Tests for the command that adds a compound skipped by `import sdf` to its wells.
"""

import shutil
import tempfile
from io import StringIO
from os.path import join

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase

from compoundlib.models import Compound
from core.models import Plate, PlateDimension, Well, WellCompound

# A structure RDKit cannot sanitize: the oxygen atom has three bonds.
BROKEN_RECORD = """Broken?Compound
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
M  END
> <CompoundName>
Broken?Compound

> <CatalogNumber>
T3608

> <WelCoordinate>
L11

> <CompoundPlateBarcode_Copy1>
TEST_A

> <CompoundPlateBarcode_Copy2>
TEST_B

> <CompoundPlateBarcode_Copy3>
TEST_MISSING

> <Vol_Copy1>
6

> <Vol_Copy2>
<24

> <Vol_Copy3>
10

$$$$
"""

# A normal record of the same library, which the command must not touch.
OTHER_RECORD = """Ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
> <CompoundName>
Ethanol

> <CatalogNumber>
T0001

> <WelCoordinate>
A3

> <CompoundPlateBarcode_Copy1>
TEST_A

> <CompoundPlateBarcode_Copy2>
TEST_B

> <CompoundPlateBarcode_Copy3>
TEST_MISSING

> <Vol_Copy1>
6

> <Vol_Copy2>
6

> <Vol_Copy3>
6

$$$$
"""

MAPPING = """compound:
  identifier: CatalogNumber
  name: CompoundName
  structure: Structure
plate:
  barcode:
    - CompoundPlateBarcode_Copy1
    - CompoundPlateBarcode_Copy2
    - CompoundPlateBarcode_Copy3
  position: WelCoordinate
  amount:
    - Vol_Copy1
    - Vol_Copy2
    - Vol_Copy3
"""

SMILES = "CC(=O)Oc1ccccc1C(=O)[O-].CC(=O)Oc1ccccc1C(=O)[O-].NC(N)=O.[Ca+2]"


class RepairLibraryWellsTest(TestCase):
    fixtures = ["plate_dimensions", "well_types"]

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.sdf_file = join(self.folder, "library.sdf")
        self.mapping_file = join(self.folder, "library_mapping.yml")
        with open(self.sdf_file, "w") as file:
            file.write(OTHER_RECORD + BROKEN_RECORD)
        with open(self.mapping_file, "w") as file:
            file.write(MAPPING)

        dimension = PlateDimension.objects.get(name="dim_384_16x24")
        self.plates = [
            Plate.objects.create(barcode=barcode, dimension=dimension)
            for barcode in ("TEST_A", "TEST_B")
        ]
        self.l11 = dimension.position("L11")

    def tearDown(self):
        shutil.rmtree(self.folder)

    def run_command(self, *extra_arguments, smiles=SMILES):
        output = StringIO()
        call_command(
            "repair_library_wells",
            "--input_file",
            self.sdf_file,
            "--mapping-file",
            self.mapping_file,
            "--identifier",
            "T3608",
            "--name",
            "Carbasalate calcium",
            "--smiles",
            smiles,
            *extra_arguments,
            stdout=output,
        )
        return output.getvalue()

    def test_the_compound_is_added_to_the_well_of_every_existing_plate(self):
        output = self.run_command()

        compound = Compound.objects.get(name="Carbasalate calcium")
        self.assertEqual(SMILES, compound.structure)
        self.assertEqual("T3608", compound.data["CatalogNumber"])
        for plate in self.plates:
            well_compound = WellCompound.objects.get(
                well__plate=plate, well__position=self.l11
            )
            self.assertEqual(compound, well_compound.compound)
            self.assertEqual(0, well_compound.amount)
        self.assertIn("Skipped: plate TEST_MISSING does not exist", output)
        # Only well L11, nothing for the normal record in A3
        self.assertEqual(2, Well.objects.count())

    def test_a_dry_run_saves_nothing(self):
        output = self.run_command("--dry-run")

        self.assertIn("Created: TEST_A L11: Carbasalate calcium", output)
        self.assertFalse(Compound.objects.exists())
        self.assertFalse(Well.objects.exists())

    def test_running_it_twice_does_not_create_duplicates(self):
        self.run_command()
        output = self.run_command()

        self.assertIn("Already there: TEST_A L11: Carbasalate calcium", output)
        self.assertEqual(1, Compound.objects.count())
        self.assertEqual(2, WellCompound.objects.count())

    def test_a_well_with_another_compound_stops_without_changes(self):
        other = Compound.objects.create(name="Something else")
        well = Well.objects.create(plate=self.plates[1], position=self.l11)
        WellCompound.objects.create(well=well, compound=other)

        with self.assertRaises(CommandError):
            self.run_command()

        self.assertFalse(Compound.objects.filter(name="Carbasalate calcium").exists())
        self.assertEqual(1, WellCompound.objects.count())

    def test_an_unreadable_smiles_is_refused(self):
        with self.assertRaises(CommandError):
            self.run_command(smiles="not a smiles")

        self.assertFalse(Compound.objects.exists())
