import json
from os.path import isfile
import re

import yaml

from .helper import sameSchema


class MappingFileSchemaError():
    def __init__(self, mapping={},
                 message="Error in mapping file schema. Should be {}"):
        self.message = message.format(mapping)
        super().__init__(self.message)


class SdfMapping():
    DEFAULT_MAPPING = {'compound': {'identifier': 'IDNUMBER', 'name': 'NAME',
        'structure': 'STRUCTURE'},
        'plate': {'barcode': 'PLATE_NUMBER1', 'position': 'POS_IN_PLATE',
            'amount': 'PLATE_AMOUNT1', }}

    def __init__(self, mappingFile: str = None):
        self.mapping = self.DEFAULT_MAPPING
        if mappingFile:
            self.load(mappingFile)

    def __str__(self):
        return json.dumps(self.mapping)

    def __repr__(self):
        return self.__str__()

    def load(self, mappingFile: str):
        if not isfile(mappingFile):
            raise FileNotFoundError(f"Cannot find mapping file {mappingFile}.")

        with open(mappingFile, 'r') as mf:
            mapping = yaml.load(mf, yaml.SafeLoader)
        if not sameSchema(self.mapping, mapping):
            raise MappingFileSchemaError(self.DEFAULT_MAPPING)

        self.mapping.update(mapping)

    @property
    def name(self) -> str:
        return self.mapping['compound']['name']

    @property
    def identifier(self) -> str:
        return self.mapping['compound']['identifier']

    @property
    def structure(self) -> str:
        return self.mapping['compound']['structure']

    @property
    def barcodes(self) -> tuple[str]:
        if isinstance(self.mapping['plate']['barcode'], str):
            return (self.mapping['plate']['barcode'],)
        return tuple(self.mapping['plate']['barcode'])

    @property
    def position(self) -> str:
        return self.mapping['plate']['position']

    @property
    def amounts(self) -> tuple[float]:
        if isinstance(self.mapping['plate']['amount'], str):
            return (self.mapping['plate']['amount'],)
        return tuple(self.mapping['plate']['amount'])


class EchoMapping():
    DEFAULT_COLUMN_HEADERS = {"source_plate_name": "Source Plate Name",
        "source_plate_barcode": "Source Plate Barcode",
        "source_plate_type": "Source Plate Type", "source_well": "Source Well",
        "source_concentration": "Source Concentration",
        "source_concentration_units": "Source Concentration Units",
        "destination_plate_name": "Destination Plate Name",
        "destination_plate_barcode": "Destination Plate Barcode",
        "destination_well": "Destination Well",
        "destination_concentration": "Destination Concentration",
        "destination_concentration_units": "Destination Concentration Units",
        "compound_name": "Compound Name", "transfer_volume": "Transfer Volume",
        "actual_volume": "Actual Volume", "transfer_status": "Transfer Status",
        "current_fluid_height": "Current Fluid Height",
        "current_fluid_volume": "Current Fluid Volume",
        "percent_dms": "% DMSO"}

    @staticmethod



    def convert_position_to_index(position: str,
                                number_of_columns: int) -> int:
        match = re.match(r'([a-zA-Z]+)(\d+)', position)
        letters = match.group(1)
        col = int(match.group(2))
        row = 0
        for char in letters:
            row += ord(char.upper()) - 65

        index = row * number_of_columns + + (col - 1)
        return index


