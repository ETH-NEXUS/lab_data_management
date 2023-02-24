import csv
import json
import os
from os.path import isfile
from typing import List, Dict
from core.models import Plate, PlateMapping
from core.mapping import Mapping, MappingList
from django.core.files.base import ContentFile

import yaml
from friendlylog import colored_logger as log

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
    def read_csv_echo_file(reader: csv.reader, headers) -> List[
        Dict[str, str]]:
        """
        Reads a CSV file and returns a list of dictionaries representing each row in the file.
        """
        table_data = []
        column_headers = list(headers.values())
        for row_index, table_row in enumerate(reader):
            if len(table_row) > 0 and table_row[0] == 'Source Plate Name':
                column_headers = table_row

            if row_index > 9:
                row_object = {}
                for cell_index, cell in enumerate(table_row):
                    row_object[column_headers[cell_index]] = cell
                table_data.append(row_object)
        return table_data

    @staticmethod
    def get_csv_echo_files(path: str, headers) -> list[
        dict[str, list[dict[str, str]] | str]]:
        """
        Recursively searches the given directory for CSV files with "transfer" in their name,
        reads each file using the read_csv_file function, and returns a list of tables
        representing the data in the CSV files.
        """
        files_data = []
        for root, dirs, files in os.walk(path, topdown=False):
            for file in files:
                if file.endswith(".csv") and 'transfer' in file:
                    subdirectory = os.path.relpath(root, path)
                    file_path = os.path.join(path, subdirectory, file)
                    with open(os.path.join(root, file), 'r') as csv_file:
                        log.debug(f"Reading file {file}")
                        csv_reader = csv.reader(csv_file, delimiter=',')
                        files_data.append({"file_path": file_path,
                                              "data": EchoMapping.read_csv_echo_file(
                                                  csv_reader, headers)})
        return files_data

    @staticmethod
    def parse_plate_data(mapping_data, file_path, headers, experiment_name):
        source_plate_barcode = mapping_data[0][headers['source_plate_barcode']]
        destination_plate_barcode = mapping_data[0][
            headers['destination_plate_barcode']]

        source_plate = Plate.objects.get(barcode=source_plate_barcode)
        destination_plate = Plate.objects.get(
            barcode=destination_plate_barcode)

        if not destination_plate:
            raise ValueError(
                f"Destination plate with barcode '{destination_plate_barcode}' "
                f"not found. Please add plates to the experiment.")
        number_of_columns = destination_plate.dimension.cols

        if destination_plate.experiment.name != experiment_name:
            raise ValueError(
                f"Destination plate with barcode '{destination_plate_barcode}' "
                f"does not belong to experiment '{experiment_name}'.")

        mapping_list = MappingList()
        for entry in mapping_data:
            if len(entry) > 0:
                log.debug(f"Mapping {entry[headers['source_well']]} to "
                          f"{entry[headers['destination_well']]} from"
                          f" {source_plate_barcode} to {destination_plate_barcode}")

                from_pos = Plate.convert_position_to_index(
                    entry[headers['source_well']], number_of_columns)
                to_pos = Plate.convert_position_to_index(
                    entry[headers['destination_well']], number_of_columns)
                mapping = Mapping(from_pos=from_pos, to_pos=to_pos,
                                  amount=float(
                                      entry[headers['actual_volume']]))
                mapping_list.add(mapping)
        mapping_success = source_plate.map(mapping_list, destination_plate)
        if mapping_success:
            file_content = open(file_path, 'rb').read()
            file_name = os.path.basename(file_path)
            plate_mapping = PlateMapping(source_plate=source_plate,
                                         target_plate=destination_plate)
            plate_mapping.mapping_file.save(file_name,
                                            ContentFile(file_content))
            log.info(f"Successfully mapped {source_plate_barcode} to "
                     f"{destination_plate_barcode}")
        else:
            log.error(f"Failed to map {source_plate_barcode} to "
                      f"{destination_plate_barcode}")
