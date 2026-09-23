"""
Reading an SDF library file with the columns its mapping file names.
Used by the `import sdf`, `fill_sdf_amounts` and `repair_library_wells` commands.
"""

from django.core.management.base import CommandError
from pandas import DataFrame
from rdkit.Chem import PandasTools

from importer.mapping import SdfMapping


def load_sdf(sdf_file: str, mapping: SdfMapping) -> DataFrame:
    """
    The records of the SDF file, one row per record and one column per
    property, e.g. {"NAME": "Aspirin", "Barcode_Copy1": "COPY_1",
    "Vol_Copy1": "24.0", "POS_IN_PLATE": "A1", "Structure": <Mol>}.

    Stops with a CommandError when a column of the mapping is missing, or when
    the number of barcode and amount columns differ.
    """
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
    check_amount_columns(mapping)
    return sdf


def check_amount_columns(mapping: SdfMapping) -> None:
    """
    The barcode and amount columns belong together by their order (the first
    amount column is the amount on the first plate copy), so there must be as
    many of each.
    """
    if len(mapping.amounts) != len(mapping.barcodes):
        raise CommandError(
            f"The mapping file names {len(mapping.barcodes)} barcode columns "
            f"but {len(mapping.amounts)} amount columns."
        )
