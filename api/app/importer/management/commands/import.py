from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from os.path import isfile
from compoundlib.models import Compound, CompoundLibrary
from core.models import Plate, Well, PlateDimension, WellCompound, WellType, Project
from platetemplate.models import PlateTemplate, PlateTemplateCategory
from importer.mapping import SdfMapping
from importer.helper import row_col_from_wells, normalize_col, normalize_row

# Excel writes an invisible BOM character at the start of a CSV file, and
# some files are not UTF-8 at all. detect_encoding finds both.
from importer.mappers.base import detect_encoding
from core.utils.plates.positions import PositionMapper
from core.models import WellDetail, PlateDetail
import numpy as np
from os.path import splitext
from pathlib import Path
from tqdm import tqdm
import argparse
import csv
from importer.helper import message

from rdkit.Chem import PandasTools
from rdkit.Chem.rdchem import Mol
from rdkit import Chem


def well_type_by_name(name: str, well_name: str) -> WellType:
    """The well type with this name, e.g. "R10"; an unknown name stops the import."""
    try:
        return WellType.by_name(name)
    except WellType.DoesNotExist:
        known_names = ", ".join(
            WellType.objects.order_by("id").values_list("name", flat=True)
        )
        raise CommandError(
            f"Unknown well type '{name}' in well {well_name}. "
            f"Known well types: {known_names}."
        )


def yes_or_no(text: str) -> bool:
    """A command line value like "yes", "no", "true" or "False" as a boolean."""
    if text.lower() in ("yes", "true", "1"):
        return True
    if text.lower() in ("no", "false", "0"):
        return False
    raise argparse.ArgumentTypeError(f"'{text}' is not yes or no")


def full_strip(text: str) -> str:
    """The text without the spaces around it."""
    return text.strip()


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "what",
            type=str,
            choices=("sdf", "template", "library_plate"),
            help="What to import",
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
            type=yes_or_no,
            help="If the plate is a control plate for an experiment: yes or no",
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

        if not isfile(sdf_file):
            raise CommandError(f"File does not exist: {sdf_file}")

        message(f"Importing SDF file {sdf_file}...", "info", room_name)

        if number_of_wells:
            number_of_rows, number_of_columns = row_col_from_wells(number_of_wells)
        if not library_name:
            library_name = splitext(Path(sdf_file).name)[0]
        library, created = CompoundLibrary.objects.update_or_create(
            name=library_name, defaults={"file_name": Path(sdf_file).name}
        )
        if created:
            message(f"Created library {library}.", "info", room_name)
        else:
            message(f"Using library {library}.", "info", room_name)

        sdf = PandasTools.LoadSDF(
            sdf_file,
            molColName=mapping.structure,
            embedProps=False,
        )
        required_columns = [mapping.name, mapping.position]
        required_columns += list(mapping.barcodes) + list(mapping.amounts)
        missing_columns = [
            column for column in required_columns if column not in sdf.columns
        ]
        if missing_columns:
            raise CommandError(
                f"These columns are not in the SDF file {sdf_file}: "
                f"{', '.join(missing_columns)}. Check the mapping file."
            )
        if len(mapping.amounts) != len(mapping.barcodes):
            raise CommandError(
                f"The mapping file names {len(mapping.barcodes)} barcode columns "
                f"but {len(mapping.amounts)} amount columns."
            )

        # Import plates
        # The plates by barcode, so the wells below do not load a plate per row
        plates = {}
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
                        defaults={"name": f"dim_{max_col*max_row}_{max_col}x{max_row}"},
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
                    plates[plate_id] = plate
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
                        mapping.name,
                    ]:
                        del data[key]

                    # The oldest one, if the same name was imported more than once
                    existing = (
                        Compound.objects.filter(name=row[mapping.name])
                        .order_by("id")
                        .first()
                    )
                    compound, created = Compound.objects.update_or_create(
                        id=existing.id if existing else None,
                        defaults={
                            "name": row[mapping.name],
                            "structure": (
                                Chem.MolToSmiles(row[mapping.structure])
                                if isinstance(row[mapping.structure], Mol)
                                else row[mapping.structure]
                            ),
                            "data": data,
                        },
                    )
                    if created:
                        __debug(f"Created compound {compound}")
                    else:
                        __debug(f"Using compound {compound}")

                    plate = plates[row[mapping_barcode]]
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
                        if isinstance(row[mapping.amounts[mapping_barcode_idx]], float)
                        or isinstance(row[mapping.amounts[mapping_barcode_idx]], int)
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
        with open(input_file, "r", encoding=detect_encoding(input_file)) as file:
            reader = csv.reader(file)
            all_rows = list(reader)
            if not any(any(cell for cell in row) for row in all_rows):
                return False, "The file is empty."
            # An editor that ends the file with a newline adds an empty last row
            while all_rows and all(cell == "" for cell in all_rows[-1]):
                all_rows.pop()
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
            columns = len(all_rows[0])
            for number, row in enumerate(before_empty_line + after_empty_line, start=1):
                if len(row) != columns:
                    return (
                        False,
                        f"Every line must have {columns} cells, like the first line, "
                        f"but line {number} has {len(row)}.",
                    )
            return True, "File format is correct."

    def __parse_library_plate_file(
        self,
        input_file,
        room_name: str = None,
    ):
        is_file_correct, message_text = self.__check_file_format(input_file)
        if not is_file_correct:
            raise CommandError(f"File format is incorrect: {message_text}")

        message("Reading plate file...", "info", room_name)
        with open(input_file, "r", encoding=detect_encoding(input_file)) as file:
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

    def library_plate(
        self,
        input_file: str,
        library_name: str,
        project_name: str,
        plate_barcode: str,
        room_name: str = None,
        is_control_plate: bool = False,
    ):
        """Imports both library plates and the control plates of a project."""
        if not isfile(input_file):
            raise CommandError(f"File does not exist: {input_file}")

        compounds, types = self.__parse_library_plate_file(input_file, room_name)
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
                raise CommandError(f"Project {project_name} does not exist.")
        else:
            raise CommandError("Please specify a library name or a project name.")

        if plate_created:
            message(
                f"Created plate {plate_barcode}.",
                "success",
                room_name,
            )

        new_compounds: list[str] = []
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
                    well_type = well_type_by_name(_type, dimension.hr_position(pos))
                    well.type = well_type
                    # The oldest one, if the same name was imported more than once
                    compound = (
                        Compound.objects.filter(name=_compound).order_by("id").first()
                    )
                    if compound is None:
                        compound = Compound.objects.create(name=_compound)
                        new_compounds.append(_compound)

                    WellCompound.objects.create(well=well, compound=compound)
                    well.save()
                pbar.update(1)

        if new_compounds:
            message(
                f"Created {len(new_compounds)} new compounds: "
                f"{', '.join(new_compounds)}",
                "success",
                room_name,
            )
        message(f"Finished processing plate {plate_barcode}.", "success", room_name)

    def template(
        self,
        input_file: str,
        category_name: str,
        template_name: str,
        room_name: str = None,
    ):
        if isfile(input_file):
            # One row of the plate per line, one well type per cell, e.g. "C<TAB>R10<TAB>P1"
            with open(input_file, "r", encoding=detect_encoding(input_file)) as file:
                message("Reading template file...", "info", room_name)
                rows = list(csv.reader(file, delimiter="\t"))
            if not rows:
                raise CommandError(f"The template file {input_file} is empty.")
            num_rows = len(rows)
            num_cols = len(rows[0])
            for row_number, row in enumerate(rows, start=1):
                if len(row) != num_cols:
                    raise CommandError(
                        f"Row {row_number} of the template file has {len(row)} cells, "
                        f"but row 1 has {num_cols}."
                    )
            well_types = [cell for row in rows for cell in row]

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
                for pos, type_name in enumerate(well_types):
                    well_name = dimension.hr_position(pos)
                    if not type_name.strip():
                        raise CommandError(
                            f"Well {well_name} has no well type in the template file."
                        )
                    well = plate.well_at(pos, create_if_not_exist=True)
                    well.type = well_type_by_name(type_name.strip(), well_name)
                    well.save()
                    pbar.update(1)
            message(
                f"Successfully imported template {template_name}.", "info", room_name
            )

        else:
            raise CommandError(f"File does not exist: {input_file}")

    def import_file(self, options: dict) -> None:
        """Imports the file of the command."""
        if not options.get("input_file"):
            raise CommandError("Please give the file to import.")
        if options.get("what") == "sdf":
            mapping = SdfMapping(options.get("mapping_file"))
            self.sdf(
                options.get("input_file"),
                mapping,
                library_name=options.get("library_name"),
                number_of_rows=options.get("number_of_rows"),
                number_of_columns=options.get("number_of_columns"),
                number_of_wells=options.get("number_of_wells"),
                debug=options.get("debug", False),
                room_name=options.get("room_name"),
            )
        elif options.get("what") == "template":
            self.template(
                options.get("input_file"),
                category_name=options.get("category_name") or "Default",
                template_name=options.get("template_name") or "Default",
                room_name=options.get("room_name"),
            )
        elif options.get("what") == "library_plate":
            if not (options.get("plate_barcode")):
                raise CommandError(
                    "Please specify plate_barcode when importing a library plate."
                )
            else:
                # check if a plate with given barcode already exists

                plate_with_given_barcode = Plate.objects.filter(
                    barcode=options.get("plate_barcode")
                )
                if plate_with_given_barcode:
                    raise CommandError(
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

    def handle(self, *args, **options):
        room_name = options.get("room_name")
        try:
            # One import is saved as a whole: after an error nothing of it is stored
            with transaction.atomic():
                self.import_file(options)
        except Exception:
            # The caller shows the error itself: the page through the Celery task,
            # the command line through the exit code
            message(
                "This import failed, nothing of it was stored. The reason:",
                "warning",
                room_name,
            )
            raise

        message("Refreshing materialized views...", "info", room_name)
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
