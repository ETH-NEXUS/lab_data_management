"""
Fills the amounts of an SDF library that is already imported, from its SDF file.

Until September 2026 the SDF import stored every amount as 0 (RDKit reads the
volumes as text, and the import only took numbers). A new import of the same
file would also change the compounds and could move plates to another library,
so this command changes only `WellCompound.amount` of the plates of one library.

Example:
    python manage.py fill_sdf_amounts --library_name 20200617_EPC_1250 \
        --input_file /data/sdf/20200617_EPC_1250.sdf \
        --mapping-file /data/sdf/20200617_EPC_1250_mapping.yml --dry-run
"""

from collections import Counter
from os.path import isfile

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from compoundlib.models import CompoundLibrary
from core.models import Plate, PlateDetail, WellCompound, WellDetail
from importer.helper import message
from importer.mapping import MappingFileSchemaError, SdfMapping
from importer.sdf_amounts import (
    amount_in_nanoliter,
    is_volume_column,
    not_a_number_warning,
    unknown_unit_warning,
)
from importer.sdf_file import is_empty_barcode, load_sdf

# How many barcodes or wells a message lists before it only counts them
LISTED_EXAMPLES = 10


def some_of(items: list) -> str:
    """The first items as text, e.g. "A, B, C and 12 more"."""
    text = ", ".join(str(item) for item in items[:LISTED_EXAMPLES])
    if len(items) > LISTED_EXAMPLES:
        text += f" and {len(items) - LISTED_EXAMPLES} more"
    return text


class Command(BaseCommand):
    help = "Fill the amounts (nL) of an imported SDF library from its SDF file."

    def add_arguments(self, parser):
        parser.add_argument("--library_name", type=str, required=True)
        parser.add_argument("--input_file", type=str, required=True)
        parser.add_argument(
            "--mapping-file",
            type=str,
            help="The mapping file of the SDF file; without it the default mapping",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Only report what would change, store nothing",
        )

    def handle(self, *args, **options):
        library = CompoundLibrary.objects.filter(
            name=options["library_name"]
        ).first()
        if library is None:
            raise CommandError(f"There is no library {options['library_name']}.")
        if not isfile(options["input_file"]):
            raise CommandError(f"File does not exist: {options['input_file']}")
        try:
            mapping = SdfMapping(options.get("mapping_file"))
        except (FileNotFoundError, MappingFileSchemaError) as error:
            raise CommandError(str(error))

        message(f"Reading SDF file {options['input_file']}...")
        sdf = load_sdf(options["input_file"], mapping)

        with transaction.atomic():
            self.fill_amounts(library, sdf, mapping)
            if options["dry_run"]:
                transaction.set_rollback(True)
                message("Dry run: nothing was stored.", "warning")
                return

        message("Refreshing materialized views...")
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)

    def fill_amounts(self, library: CompoundLibrary, sdf, mapping: SdfMapping) -> None:
        plates = {
            plate.barcode: plate
            for plate in Plate.objects.filter(library=library).select_related(
                "dimension"
            )
        }
        # The compounds of the library wells by plate barcode, well position and
        # compound name, e.g. {("EPC_1250_1A", 0, "Aspirin"): <WellCompound>}
        well_compounds = {
            (
                well_compound.well.plate.barcode,
                well_compound.well.position,
                well_compound.compound.name,
            ): well_compound
            for well_compound in WellCompound.objects.filter(
                well__plate__library=library
            ).select_related("well__plate", "compound")
        }

        changed = []
        unchanged = 0
        empty_barcodes = 0
        plates_of_other_libraries = set()
        plates_in_file = set()
        # e.g. ["EPC_1250_1A A1 Aspirin"]
        missing_well_compounds = []

        for barcode_column, amount_column in zip(mapping.barcodes, mapping.amounts):
            if not is_volume_column(amount_column):
                message(unknown_unit_warning(amount_column), "warning")
                continue
            # The values that are not a number, e.g. {"<24": 2486}
            not_a_number: Counter[str] = Counter()

            for _, row in sdf.iterrows():
                barcode = row[barcode_column]
                if is_empty_barcode(barcode):
                    empty_barcodes += 1
                    continue
                plate = plates.get(barcode)
                if plate is None:
                    plates_of_other_libraries.add(barcode)
                    continue
                plates_in_file.add(barcode)

                position = plate.dimension.position(row[mapping.position])
                well_compound = well_compounds.get(
                    (barcode, position, row[mapping.name])
                )
                if well_compound is None:
                    missing_well_compounds.append(
                        f"{barcode} {row[mapping.position]} {row[mapping.name]}"
                    )
                    continue

                amount = amount_in_nanoliter(row[amount_column])
                if amount is None:
                    not_a_number[str(row[amount_column])] += 1
                    amount = 0.0
                if well_compound.amount == amount:
                    unchanged += 1
                else:
                    well_compound.amount = amount
                    changed.append(well_compound)

            if not_a_number:
                message(not_a_number_warning(amount_column, not_a_number), "warning")

        WellCompound.objects.bulk_update(changed, ["amount"], batch_size=1000)

        self.report(
            library,
            changed=len(changed),
            unchanged=unchanged,
            empty_barcodes=empty_barcodes,
            plates_of_other_libraries=sorted(plates_of_other_libraries),
            plates_not_in_file=sorted(set(plates) - plates_in_file),
            missing_well_compounds=missing_well_compounds,
        )

    def report(
        self,
        library: CompoundLibrary,
        changed: int,
        unchanged: int,
        empty_barcodes: int,
        plates_of_other_libraries: list,
        plates_not_in_file: list,
        missing_well_compounds: list,
    ) -> None:
        message(
            f"Library {library.name}: {changed} amounts changed, "
            f"{unchanged} already had the amount of the file.",
            "success",
        )
        if empty_barcodes:
            message(
                f"{empty_barcodes} plate copies of the file have no barcode "
                f"and were skipped.",
            )
        if plates_of_other_libraries:
            message(
                f"{len(plates_of_other_libraries)} plates of the file are not "
                f"in library {library.name} and were skipped: "
                f"{some_of(plates_of_other_libraries)}.",
                "warning",
            )
        if plates_not_in_file:
            message(
                f"{len(plates_not_in_file)} plates of library {library.name} are "
                f"not in the file (with this mapping file): "
                f"{some_of(plates_not_in_file)}.",
                "warning",
            )
        if missing_well_compounds:
            message(
                f"{len(missing_well_compounds)} compounds of the file are not in "
                f"their well in LDM and were skipped: "
                f"{some_of(missing_well_compounds)}.",
                "warning",
            )
