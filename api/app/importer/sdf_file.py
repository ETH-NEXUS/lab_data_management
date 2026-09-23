"""
Reading an SDF library file with the columns its mapping file names.
Used by the `import sdf` and the `fill_sdf_amounts` commands.
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
    the number of barcode and amount columns differ (they belong together by
    their order: the first amount column is the amount on the first plate copy).
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
    if len(mapping.amounts) != len(mapping.barcodes):
        raise CommandError(
            f"The mapping file names {len(mapping.barcodes)} barcode columns "
            f"but {len(mapping.amounts)} amount columns."
        )
    return sdf
