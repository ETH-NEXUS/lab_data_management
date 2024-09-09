from django.core.management.base import BaseCommand
import traceback
from os.path import isfile
from compoundlib.models import Compound, CompoundLibrary
from core.models import Plate, Well, PlateDimension, WellCompound, WellType, Project
from platetemplate.models import PlateTemplate, PlateTemplateCategory
from importer.mapping import SdfMapping
from importer.helper import row_col_from_wells, normalize_col, normalize_row
from core.mapping import PositionMapper
from core.models import WellDetail, PlateDetail
import numpy as np
from os.path import splitext
from pathlib import Path
from tqdm import tqdm
import csv
from importer.helper import message

from rdkit.Chem import PandasTools
from rdkit.Chem.rdchem import Mol
from rdkit import Chem


def full_strip(s: str):
    return s.strip().lstrip()


class Command(BaseCommand):
    def add_arguments(self, parser):

        parser.add_argument(
            "--input_file",
            "-i",
            type=str,
            required=True,
            help="The input file",
        )
        parser.add_argument(
            "--mapping-file",
            "-m",
            type=str,
            help="The mapping file for the sdf columns, otherwise default mapping is used",
        )

    def compound_data(
        self,
    ):
        pass

    def handle(self, *args, **options):
        try:
            imported = False

            if imported:
                message(
                    "Refreshing materialized views...", "info", options.get("room_name")
                )
                PlateDetail.refresh(concurrently=True)
                WellDetail.refresh(concurrently=True)

        except Exception as ex:
            message(ex, "error", options.get("room_name"))
            traceback.print_exc()
