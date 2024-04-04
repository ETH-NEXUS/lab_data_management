from django.core.management.base import BaseCommand
import csv
import pandas as pd
from tqdm import tqdm
from helpers.logger import logger

from compoundlib.models import Compound, CompoundLibrary
from core.models import (
    Plate,
    PlateDimension,
    Well,
    WellCompound,
    PlateDetail,
    WellDetail,
)

PLACEHOLDER_STRUCTURE = "COC(=O)C(Cc1ccccc1)NC(=O)C(CC(C)C)NC(=O)OC(C)(C)C"
DIMENSION = PlateDimension.objects.get(name="dim_1536_32x48")


def read_file(file_path):
    with open(file_path, "r") as file:
        logger.info(f"Reading file {file_path}")
        return list(csv.DictReader(file, delimiter="\t"))


def read_xlsx(file_path):
    logger.info(f"Reading xlsx file {file_path}")
    df = pd.read_excel(file_path)
    return df.to_dict(orient="records")


class Command(BaseCommand):
    help = "Import fixed plate data from a file into the database."

    def add_arguments(self, parser):
        parser.add_argument(
            "what",
            type=str,
            help="What to do",
            choices=["change_plate", "import_xlsx"],
        )
        parser.add_argument(
            "-i", "--input_file", required=True, help="Path to the input file."
        )
        parser.add_argument(
            "-l", "--library_name", required=True, help="Name of the compound library."
        )

    def handle(self, *args, **options):
        what = options["what"]
        library = self.get_or_create_library(options["library_name"])
        input_file = options["input_file"]
        if what == "change_plate":
            rows = read_file(input_file)
            self.process_rows(rows, library)
        elif what == "import_xlsx":
            self.import_xlsx(input_file, library)

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    def import_xlsx(self, file_path, library):
        data = read_xlsx(file_path)
        plate_cache = {}
        well_compound_batch = []

        with tqdm(total=len(data)) as pbar:
            for item in data:
                compound = self.get_or_create_compound(
                    item, use_placeholder=True, name=item["Compound_Id"]
                )

                plate_barcode = item["Labware_Barcode"]
                if plate_barcode not in plate_cache:
                    plate, _ = Plate.objects.get_or_create(
                        barcode=plate_barcode,
                        library=library,
                        defaults={"dimension": DIMENSION},
                    )
                    plate_cache[plate_barcode] = plate
                else:
                    plate = plate_cache[plate_barcode]

                well = self.get_or_create_well(item["Well"], plate)
                well_compound = WellCompound(compound=compound, well=well)
                well_compound_batch.append(well_compound)
                pbar.update(1)

        WellCompound.objects.bulk_create(well_compound_batch, ignore_conflicts=True)

    def get_or_create_library(self, library_name):
        library, _ = CompoundLibrary.objects.get_or_create(name=library_name)
        logger.info(f"Library {library.name} found.")
        return library

    def process_rows(self, rows, library):
        unique_barcodes = self.extract_unique_barcodes(rows)
        logger.info(f"Unique barcodes: {unique_barcodes}")
        new_plates = self.create_plates(unique_barcodes, library)
        self.create_well_compounds(rows, new_plates)

    def extract_unique_barcodes(self, rows):
        barcodes = []
        for row in rows:
            for key, value in row.items():
                if key.startswith("CompoundPlateBarcode"):
                    barcodes.append(value)
        return list(set(barcodes))

    def create_plates(self, barcodes, library):
        dimension_384 = PlateDimension.objects.get(name="dim_384_16x24")
        new_plates = []
        for barcode in barcodes:
            new_barcode = f"{barcode}_fixed"
            plate, _ = Plate.objects.get_or_create(
                barcode=new_barcode, library=library, dimension=dimension_384
            )
            plate.save()
            if _:
                logger.info(f"Plate {plate.barcode} created.")
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

    def get_or_create_compound(self, row, use_placeholder=False, name=None):
        if use_placeholder:
            compound, _ = Compound.objects.get_or_create(
                name=name, structure=PLACEHOLDER_STRUCTURE
            )
            compound.save()
            return compound
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
