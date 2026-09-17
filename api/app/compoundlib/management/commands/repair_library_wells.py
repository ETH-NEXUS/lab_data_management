"""
Adds one compound that `import sdf` skipped to the wells of its library plates.

`import sdf` drops a whole SDF record when RDKit cannot read its structure. Then
the compound and its well are missing in every plate copy of the record.
Example: "Carbasalate?Calcium" in 20211119_NEXUS_FDA_approved.sdf has an oxygen
atom with three bonds, so well L11 is empty in Drug08_A ... Drug08_I.

This command reads that one record from the SDF file without checking the
structure, takes the well and the plate barcodes from it, and creates the
compound, the well and the well compound. The corrected structure and the name
are given on the command line.

Example arguments (run with --dry-run first, then without it):

    python manage.py repair_library_wells \
        --input_file /data/sdf/20211119_NEXUS_FDA_approved.sdf \
        --mapping-file /data/sdf/20211119_NEXUS_FDA_approved_mapping.yml \
        --identifier T3608 \
        --name "Carbasalate calcium" \
        --smiles "CC(=O)Oc1ccccc1C(=O)[O-].CC(=O)Oc1ccccc1C(=O)[O-].NC(N)=O.[Ca+2]" \
        --dry-run

--identifier is the value of the mapping's compound identifier column
(CatalogNumber in the FDA approved mapping file).
"""

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from rdkit import Chem, RDLogger

from compoundlib.models import Compound
from core.models import (
    ExperimentDetail,
    Plate,
    PlateDetail,
    Well,
    WellCompound,
    WellDetail,
)
from importer.mapping import SdfMapping

# Same as the other wells of the library: `import sdf` stores 0, because RDKit
# reads the volumes of an SDF file as text.
WELL_COMPOUND_AMOUNT = 0


def read_sdf_record(sdf_file: str, identifier_column: str, identifier: str) -> dict:
    """
    Returns the fields of the one record whose identifier column has the value.
    The structure is not checked, so records that `import sdf` skipped are found.

    Returned data example:
    {"CatalogNumber": "T3608", "WelCoordinate": "L11",
     "CompoundPlateBarcode_Copy1": "Drug08_A", "Vol_Copy1": "6", ...}
    """
    # RDKit prints a warning for every structure it cannot sanitize
    RDLogger.DisableLog("rdApp.*")
    try:
        records = []
        for molecule in Chem.SDMolSupplier(sdf_file, sanitize=False):
            if molecule is None:
                continue
            fields = {name: molecule.GetProp(name) for name in molecule.GetPropNames()}
            if fields.get(identifier_column) == identifier:
                records.append(fields)
    finally:
        RDLogger.EnableLog("rdApp.*")

    if len(records) != 1:
        raise CommandError(
            f"Expected one record with {identifier_column} = {identifier} "
            f"in {sdf_file}, found {len(records)}."
        )
    return records[0]


class Command(BaseCommand):
    help = "Add a compound that the SDF import skipped to its library plate wells."

    def add_arguments(self, parser):
        parser.add_argument(
            "--input_file", "-i", required=True, help="The SDF file of the library"
        )
        parser.add_argument(
            "--mapping-file",
            "-m",
            required=True,
            help="The mapping file that was used for `import sdf`",
        )
        parser.add_argument(
            "--identifier",
            required=True,
            help="Value of the compound identifier column, e.g. T3608",
        )
        parser.add_argument("--name", required=True, help="Name of the compound")
        parser.add_argument(
            "--smiles", required=True, help="The corrected structure as SMILES"
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Show what would be created, without saving anything",
        )

    def handle(self, *args, **options):
        molecule = Chem.MolFromSmiles(options["smiles"])
        if molecule is None:
            raise CommandError(f"RDKit cannot read the SMILES {options['smiles']}.")

        mapping = SdfMapping(options["mapping_file"])
        fields = read_sdf_record(
            options["input_file"], mapping.identifier, options["identifier"]
        )
        well_coordinate = fields[mapping.position]

        with transaction.atomic():
            compound, created = Compound.objects.get_or_create(
                name=options["name"],
                defaults={
                    "structure": Chem.MolToSmiles(molecule),
                    # Like `import sdf`: the other fields of the record
                    "data": {
                        key: value
                        for key, value in fields.items()
                        if key not in (mapping.name, mapping.position)
                    },
                },
            )
            self.report("Created" if created else "Using", f"compound {compound}")

            for barcode_column in mapping.barcodes:
                barcode = fields.get(barcode_column)
                if barcode:
                    self.add_to_plate(barcode, well_coordinate, compound)

            if options["dry_run"]:
                # Everything above ran for real, but nothing is kept
                transaction.set_rollback(True)
                self.report("Dry run", "nothing was saved")
                return

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)
        self.report("Done", "materialized views refreshed")

    def add_to_plate(self, barcode: str, well_coordinate: str, compound: Compound):
        plate = Plate.objects.filter(barcode=barcode).first()
        if plate is None:
            self.report("Skipped", f"plate {barcode} does not exist")
            return

        position = plate.dimension.position(well_coordinate)
        well, created = Well.objects.get_or_create(plate=plate, position=position)

        other_compounds = well.well_compounds.exclude(compound=compound)
        if other_compounds.exists():
            # Stops the whole command: the transaction undoes all changes
            names = ", ".join(str(item.compound) for item in other_compounds)
            raise CommandError(
                f"Well {well_coordinate} of plate {barcode} already contains "
                f"{names}. Nothing was changed."
            )

        _, well_compound_created = WellCompound.objects.get_or_create(
            well=well,
            compound=compound,
            defaults={"amount": WELL_COMPOUND_AMOUNT},
        )
        self.report(
            "Created" if well_compound_created else "Already there",
            f"{barcode} {well_coordinate}: {compound}",
        )

    def report(self, action: str, text: str):
        self.stdout.write(f"{action}: {text}")
