from django.core.management.base import BaseCommand
import traceback
from friendlylog import colored_logger as log
from os.path import isfile
from compoundlib.models import Compound, CompoundLibrary
from core.models import Plate, Well, PlateDimension
from importer.mapping import SdfMapping
from core.mapping import PositionMapper
from pprint import pprint as pp
from copy import deepcopy
import numpy as np
from os.path import splitext
from pathlib import Path

from rdkit.Chem import PandasTools
from rdkit.Chem.rdchem import Mol
from rdkit import Chem


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('action', type=str, help='The action to execute')
        parser.add_argument('--input-file', '-i', type=str, required=True, help='The input file')
        parser.add_argument('--mapping-file', '-m', type=str, help='The mapping file for the sdf columns, otherwise default mapping is used')
        parser.add_argument('--library-name', '-n', type=str, help='The name of the library, otherwise the filename is the library name')

    def sdf(self, sdf_file, mapping: SdfMapping, library_name=None):
        log.info(f"Importing SDF file {sdf_file}...")
        if isfile(sdf_file):
            if not library_name:
                library_name = splitext(Path(sdf_file).name)[0]
            library, created = CompoundLibrary.objects.update_or_create(
                name=library_name,
                defaults={
                    'file_name': Path(sdf_file).name
                }
            )
            if created:
                log.debug(f"Created library {library}.")
            else:
                log.debug(f"Using library {library}.")

            sdf = PandasTools.LoadSDF(sdf_file, molColName='Structure', smilesName='NAME', embedProps=False,
                                      includeFingerprints=True)

            # Import plates
            for mapping_barcode_idx, mapping_barcode in enumerate(mapping.barcodes):
                for plate_id in sdf[mapping_barcode].unique():
                    # Determinate Plate Dimension
                    max_row = 0
                    max_col = 0
                    for position in sdf.loc[sdf[mapping_barcode] == plate_id][mapping.position]:
                        row, col = PositionMapper.map(position)
                        max_row = max(max_row, row)
                        max_col = max(max_col, col)

                    plateDimension, created = PlateDimension.objects.get_or_create(
                        rows=max_row,
                        cols=max_col,
                        defaults={
                            'name': f"dim_{max_col}x{max_row}"
                        })
                    if created:
                        log.debug(
                            f"Created plate dimension {plateDimension}.")
                    else:
                        log.debug(f"Using plate dimension {plateDimension}.")

                    plate, created = Plate.objects.update_or_create(
                        barcode=plate_id,
                        defaults={
                            'dimension': plateDimension,
                            'library': library
                        }
                    )
                    if created:
                        log.debug(f"Created plate {plate.barcode}.")
                    else:
                        log.debug(f"Using plate {plate.barcode}.")

                # Import Compounds and Wells
                for _, row in sdf.iterrows():

                    data = row.replace({np.nan: None}).to_dict()
                    for key in [mapping.structure, mapping_barcode, mapping.amounts[mapping_barcode_idx], mapping.position, mapping.identifier]:
                        del data[key]

                    compound, created = Compound.objects.update_or_create(
                        identifier=row[mapping.identifier],
                        defaults={
                            'library': library,
                            'structure': Chem.MolToSmiles(row[mapping.structure]) if isinstance(row[mapping.structure], Mol) else row[mapping.structure],
                            'smile': row[mapping.smile],
                            'data': data
                        }
                    )
                    if created:
                        log.debug(
                            f"Created compound {compound}")
                    else:
                        log.debug(
                            f"Using compound {compound}")

                    plate = Plate.objects.get(barcode=row[mapping_barcode])
                    well, created = Well.objects.update_or_create(
                        plate=plate,
                        position=plate.dimension.position(
                            row[mapping.position]
                        ),
                        defaults={
                            'amount': row[mapping.amounts[mapping_barcode_idx]],
                            'compound': compound
                        }
                    )
                    if created:
                        log.debug(
                            f"Created well {well.plate}: {well.hr_position} ({row[mapping.position]})")
                    else:
                        log.debug(
                            f"Using well {well.plate}: {well.hr_position} ({row[mapping.position]})")

    def handle(self, *args, **options):
        try:
            if options.get('action') == 'sdf':
                mapping = SdfMapping(options.get('mapping_file'))
                self.sdf(options.get('input_file'), mapping, library_name=options.get('library_name'))

        except Exception as ex:
            log.error(ex)
            traceback.print_exc()
