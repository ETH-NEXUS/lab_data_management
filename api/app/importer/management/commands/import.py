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
            "what",
            type=str,
            help="What to import: sdf | template | library_plate",
        )
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
        parser.add_argument(
            "--library_name",
            "-l",
            type=str,
            help="The name of the library, otherwise the filename is the library name",
            default=None,
        )
        parser.add_argument(
            "--number-of-rows",
            "-r",
            type=int,
            help="The number of rows of each plate",
        )
        parser.add_argument(
            "--number-of-columns",
            "-c",
            type=int,
            help="The number of columns of each plate",
        )
        parser.add_argument(
            "--number-of-wells",
            "-n",
            type=int,
            choices=[96, 384, 1536],
            help="The number of wells each plate. This setting overrides --number-of-rows and --number-of-columns",
        )
        parser.add_argument(
            "--category_name",
            "--cat",
            type=str,
            help="The name of the template category, defaults to 'Default'",
            default="Default",
        )
        parser.add_argument(
            "--template_name",
            "-t",
            type=str,
            help="The name of the template, defaults to 'Default'",
            default="Default",
        )
        parser.add_argument(
            "--debug", action="store_true", help="Outputs debug messages"
        )
        parser.add_argument(
            "--room_name",
            "-o",
            help="Unique room name for long polling.",
        )

        parser.add_argument(
            "--plate_barcode",
            help="The barcode of the new library plate",
        )
        parser.add_argument(
            "--project_name",
            help="If you are importing a control plate, it should be associated with a project",
            default=None,
        )
        parser.add_argument(
            "--is_control_plate",
            help="If the plate is a control plate for an experiment",
            default=False,
        )

    def sdf(
        self,
        sdf_file: str,
        mapping: SdfMapping,
        library_name: str = None,
        number_of_rows: int = None,
        number_of_columns: int = None,
        number_of_wells: int = None,
        debug: bool = False,
        room_name: str = None,
    ):
        def __debug(msg):
            if debug:
                message(msg, "debug", room_name)

        message(f"Importing SDF file {sdf_file}...", "info", room_name)

        if number_of_wells:
            number_of_rows, number_of_columns = row_col_from_wells(number_of_wells)
        if isfile(sdf_file):
            if not library_name:
                library_name = splitext(Path(sdf_file).name)[0]
            library, created = CompoundLibrary.objects.update_or_create(
                name=library_name, defaults={"file_name": Path(sdf_file).name}
            )
            if created:
                message(f"Created library {library}.", "info", room_name)
                __debug(f"Created library {library}.")
            else:
                message(f"Created library {library}.", "info", room_name)
                __debug(f"Using library {library}.")

            sdf = PandasTools.LoadSDF(
                sdf_file,
                molColName=mapping.structure,
                embedProps=False,
                includeFingerprints=True,
            )

            # Import plates
            for mapping_barcode_idx, mapping_barcode in enumerate(mapping.barcodes):
                message(
                    f"Processing plates for barcode column {mapping_barcode}...",
                    "info",
                    room_name,
                )

                with tqdm(
                    desc="Processing plates",
                    unit="plates",
                    total=len(sdf[mapping_barcode].unique()),
                ) as pbar:
                    for plate_id in sdf[mapping_barcode].unique():
                        # Determinate Plate Dimension
                        if number_of_columns and number_of_rows:
                            max_row = number_of_rows
                            max_col = number_of_columns
                        else:
                            max_row = 0
                            max_col = 0
                            for position in sdf.loc[sdf[mapping_barcode] == plate_id][
                                mapping.position
                            ]:
                                row, col = PositionMapper.map(position)
                                max_row = max(max_row, row)
                                max_col = max(max_col, col)
                            max_row = normalize_row(max_row)
                            max_col = normalize_col(max_col)

                        (
                            plateDimension,
                            created,
                        ) = PlateDimension.objects.get_or_create(
                            rows=max_row,
                            cols=max_col,
                            defaults={
                                "name": f"dim_{max_col*max_row}_{max_col}x{max_row}"
                            },
                        )
                        if created:
                            __debug(f"Created plate dimension {plateDimension}.")
                        else:
                            __debug(f"Using plate dimension {plateDimension}.")

                        plate, created = Plate.objects.update_or_create(
                            barcode=plate_id,
                            defaults={
                                "dimension": plateDimension,
                                "library": library,
                            },
                        )
                        if created:
                            __debug(f"Created plate {plate.barcode}.")
                        else:
                            __debug(f"Using plate {plate.barcode}.")
                        pbar.update(1)

                # Import Compounds and Wells
                with tqdm(
                    desc="Processing wells", unit="wells", total=len(sdf.index)
                ) as wbar:
                    for _, row in sdf.iterrows():
                        data = row.replace({np.nan: None}).to_dict()
                        for key in [
                            mapping.structure,
                            mapping_barcode,
                            mapping.amounts[mapping_barcode_idx],
                            mapping.position,
                            mapping.identifier,
                            mapping.name,
                        ]:
                            del data[key]

                        compound, created = Compound.objects.update_or_create(
                            identifier=row[mapping.identifier],
                            defaults={
                                "name": row[mapping.name],
                                "structure": Chem.MolToSmiles(row[mapping.structure])
                                if isinstance(row[mapping.structure], Mol)
                                else row[mapping.structure],
                                "data": data,
                            },
                        )
                        if created:
                            __debug(f"Created compound {compound}")
                        else:
                            __debug(f"Using compound {compound}")
                        compound.library = library
                        compound.save()

                        plate = Plate.objects.get(barcode=row[mapping_barcode])
                        well, created = Well.objects.update_or_create(
                            plate=plate,
                            position=plate.dimension.position(row[mapping.position]),
                        )
                        if created:
                            __debug(
                                f"Created well {well.plate}: {well.hr_position} ({row[mapping.position]})"
                            )
                        else:
                            __debug(
                                f"Using well {well.plate}: {well.hr_position} ({row[mapping.position]})"
                            )
                        amount = (
                            row[mapping.amounts[mapping_barcode_idx]]
                            if isinstance(
                                row[mapping.amounts[mapping_barcode_idx]], float
                            )
                            or isinstance(
                                row[mapping.amounts[mapping_barcode_idx]], int
                            )
                            else 0
                        )
                        (
                            well_compound,
                            created,
                        ) = WellCompound.objects.update_or_create(
                            well=well,
                            compound=compound,
                            defaults={
                                "amount": amount,
                            },
                        )
                        if created:
                            __debug(
                                f"Created well_compound {well_compound.well} -> {well_compound.compound}"
                            )
                        else:
                            __debug(
                                f"Using well_compound {well_compound.well} -> {well_compound.compound}"
                            )
                        wbar.update(1)

    def __check_file_format(self, input_file: str):
        with open(input_file, "r") as file:
            reader = csv.reader(file)
            all_rows = list(reader)
            empty_rows = [row for row in all_rows if all(x == "" for x in row)]
            if len(empty_rows) != 1:
                return False, "The file should contain exactly one empty line."
            split_index = all_rows.index(empty_rows[0])
            before_empty_line = all_rows[:split_index]
            after_empty_line = all_rows[split_index + 1 :]
            if len(before_empty_line) != len(after_empty_line):
                return (
                    False,
                    "The number of lines before and after the empty line should be equal.",
                )
            return True, "File format is correct."

    def __parse_library_plate_file(
        self,
        input_file,
        room_name: str = None,
    ):
        is_file_correct, message_text = self.__check_file_format(input_file)
        if not is_file_correct:
            message(message_text, "error", room_name)
            return None, None

        if isfile(input_file):
            message("Reading plate file...", "info", room_name)
            with open(input_file, "r") as file:
                reader = csv.reader(file)
                matrix1 = []
                matrix2 = []
                current_matrix = matrix1
                for row in reader:
                    if all(x == "" for x in row):
                        row = None
                    if not row:
                        current_matrix = matrix2
                    else:
                        current_matrix.append(row)
                return matrix1, matrix2

    """
    This function imports both library plates and project control  plates.
    """

    def library_plate(
        self,
        input_file: str,
        library_name: str,
        project_name: str,
        plate_barcode: str,
        room_name: str = None,
        is_control_plate: bool = False,
    ):
        if isfile(input_file):
            compounds, types = self.__parse_library_plate_file(input_file)
            if compounds is None or types is None:
                message("File format is incorrect.", "error", room_name)
                return
            num_cols = len(compounds[0])
            num_rows = len(compounds)
            dimension, _ = PlateDimension.objects.get_or_create(
                rows=num_rows,
                cols=num_cols,
                defaults={"name": f"dim_{num_cols*num_rows}_{num_cols}x{num_rows}"},
            )

            compounds = [item for sublist in compounds for item in sublist]
            types = [item for sublist in types for item in sublist]
            plate_created = False
            plate = None
            library = None
            if library_name:
                library, library_created = CompoundLibrary.objects.update_or_create(
                    name=library_name
                )
                if library_created:
                    message(f"Created library {library_name}.", "success", room_name)
                plate, plate_created = Plate.objects.update_or_create(
                    barcode=plate_barcode,
                    dimension=dimension,
                    library=library,
                    is_control_plate=is_control_plate,
                )
            elif project_name:
                try:
                    project = Project.objects.get(name=project_name)
                    plate, plate_created = Plate.objects.update_or_create(
                        barcode=plate_barcode,
                        dimension=dimension,
                        project=project,
                        is_control_plate=is_control_plate,
                    )
                except Project.DoesNotExist:
                    message(
                        f"Project {project_name} does not exist.", "error", room_name
                    )

            if plate_created:
                message(
                    f"Created plate {plate_barcode}.",
                    "success",
                    room_name,
                )

            with tqdm(
                desc="Processing wells",
                unit="wells",
                total=len(compounds),
            ) as pbar:
                for pos, content in enumerate(compounds):
                    _type = full_strip(types[pos])
                    _compound = full_strip(content)
                    if _type != "null" and _compound != "null":
                        well = plate.well_at(pos, create_if_not_exist=True)
                        well_type = WellType.by_name(_type)
                        well.type = well_type
                        compound, created = Compound.objects.get_or_create(
                            name=_compound,
                            library=library,
                            identifier=f"{_compound}_{library.name if library else project_name}",
                        )
                        if created:
                            message(
                                f"Created compound {_compound}.", "success", room_name
                            )
                        WellCompound.objects.create(well=well, compound=compound)
                        well.save()
                        pbar.update(1)
                print("**********************************")
            message(f"Finished processing plate {plate_barcode}.", "success", room_name)

    def template(
        self,
        input_file: str,
        category_name: str,
        template_name: str,
        room_name: str = None,
    ):
        if isfile(input_file):
            well_types = []
            num_rows = 0
            num_cols = 0
            with open(input_file, "r") as file:
                message("Reading template file...", "info", room_name)
                reader = csv.reader(file, delimiter="\t")
                for row in reader:
                    well_types += row
                    num_rows += 1
                num_cols = len(row)

            dimension, _ = PlateDimension.objects.get_or_create(
                rows=num_rows,
                cols=num_cols,
                defaults={"name": f"dim_{num_cols*num_rows}_{num_cols}x{num_rows}"},
            )

            category, _ = PlateTemplateCategory.objects.get_or_create(
                name=category_name
            )

            template, _ = PlateTemplate.objects.get_or_create(
                name=template_name, category=category
            )

            plate, _ = Plate.objects.get_or_create(
                barcode=f"__TEMPL__{category_name}_{template_name}",
                dimension=dimension,
                template=template,
            )

            with tqdm(
                desc="Processing wells",
                unit="wells",
                total=len(well_types),
            ) as pbar:
                for pos, type in enumerate(well_types):
                    well = plate.well_at(pos, create_if_not_exist=True)
                    well_type = WellType.by_name(type[0])
                    well.type = well_type
                    well.save()
                    pbar.update(1)
            message(
                f"Successfully imported template {template_name}.", "info", room_name
            )

        else:
            message(f"File does not exist: {input_file}", "error", room_name)

    def handle(self, *args, **options):

        try:
            imported = False
            if options.get("what") == "sdf":
                mapping = SdfMapping(options.get("mapping_file"))
                self.sdf(
                    options.get("input_file"),
                    mapping,
                    library_name=options.get("library_name"),
                    number_of_rows=options.get("number_of_rows"),
                    number_of_columns=options.get("number_of_columns"),
                    number_of_wells=options.get("number_of_wells"),
                    room_name=options.get("room_name"),
                )
                imported = True
            elif options.get("what") == "template":
                self.template(
                    options.get("input_file"),
                    category_name=options.get("category_name"),
                    template_name=options.get("template_name"),
                    room_name=options.get("room_name"),
                )
                imported = True
            elif options.get("what") == "library_plate":
                if not (options.get("plate_barcode")):
                    message(
                        "ERROR: Please specify plate_barcode.",
                        "error",
                        options.get("room_name"),
                    )
                    raise ValueError(
                        "Please specify plate_barcode when importing a library plate."
                    )
                else:
                    # check if a plate with given barcode already exists

                    plate_with_given_barcode = Plate.objects.filter(
                        barcode=options.get("plate_barcode")
                    )
                    if plate_with_given_barcode:
                        message(
                            f"ERROR: Plate with barcode {options.get('plate_barcode')} already exists.",
                            "error",
                            options.get("room_name"),
                        )
                        raise ValueError(
                            f"Plate with barcode {options.get('plate_barcode')} already exists."
                        )

                    self.library_plate(
                        options.get("input_file"),
                        library_name=options.get("library_name"),
                        plate_barcode=options.get("plate_barcode"),
                        room_name=options.get("room_name"),
                        is_control_plate=options.get("is_control_plate"),
                        project_name=options.get("project_name"),
                    )
                    imported = True

            if imported:
                message(
                    "Refreshing materialized views...", "info", options.get("room_name")
                )
                PlateDetail.refresh(concurrently=True)
                WellDetail.refresh(concurrently=True)

        except Exception as ex:
            message(ex, "error", options.get("room_name"))
            traceback.print_exc()
