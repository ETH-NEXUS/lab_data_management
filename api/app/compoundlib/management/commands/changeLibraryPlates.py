from django.core.management.base import BaseCommand
import csv

from compoundlib.models import Compound, CompoundLibrary
from core.models import (
    Plate,
    PlateDimension,
    Well,
    WellCompound,
    PlateDetail,
    WellDetail,
)


def read_file(file_path):
    with open(file_path, "r") as file:
        print(f"Reading file {file_path}")
        return list(csv.DictReader(file, delimiter="\t"))


class Command(BaseCommand):
    help = "Import fixed plate data from a file into the database."

    def add_arguments(self, parser):
        parser.add_argument(
            "-i", "--input_file", required=True, help="Path to the input file."
        )
        parser.add_argument(
            "-l", "--library_name", required=True, help="Name of the compound library."
        )

    def handle(self, *args, **options):
        library = self.get_or_create_library(options["library_name"])
        rows = read_file(options["input_file"])
        self.process_rows(rows, library)
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    def get_or_create_library(self, library_name):
        library, _ = CompoundLibrary.objects.get_or_create(name=library_name)
        print(f"Library {library.name} found.")
        return library

    def process_rows(self, rows, library):
        unique_barcodes = self.extract_unique_barcodes(rows)
        new_plates = self.create_plates(unique_barcodes, library)
        self.create_well_compounds(rows, new_plates)

    def extract_unique_barcodes(self, rows):
        return list(
            {
                row["CompoundPlateBarcode"]
                for row in rows
                if "CompoundPlateBarcode" in row
            }
        )

    def create_plates(self, barcodes, library):
        dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        new_plates = []
        for barcode in barcodes:
            new_barcode = f"{barcode}_fixed"
            plate, _ = Plate.objects.get_or_create(
                barcode=new_barcode, library=library, dimension=dimension_384
            )
            plate.save()
            print(f"Plate {plate.barcode} created.")
            new_plates.append(plate)
        return new_plates

    def create_well_compounds(self, rows, plates):
        for row in rows:
            compound = self.get_or_create_compound(row)
            for plate in plates:
                well = self.get_or_create_well(row["WellCoordinate"], plate)
                well_compound, created = WellCompound.objects.get_or_create(
                    compound=compound, well=well
                )
                well_compound.save()

    def get_or_create_compound(self, row):
        structure = row["Smiles"] if row["Smiles"] else None
        compound, _ = Compound.objects.get_or_create(
            name=row["MoleculeName"], structure=structure
        )

        compound.save()
        return compound

    def get_or_create_well(self, position_raw, plate):
        position = plate.dimension.position(position_raw)
        well, _ = Well.objects.get_or_create(plate=plate, position=position)
        well.save()
        return well
