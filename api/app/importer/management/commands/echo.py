import csv
import os
import re
import traceback
from typing import List, Dict
from core.mapping import Mapping, MappingList
from core.models import Plate, PlateMapping
from django.core.files import File
from django.core.management.base import BaseCommand
from friendlylog import colored_logger as log
import yaml
from django.core.files.base import ContentFile

# add plates to the experiment, put the echo files to /data/echo and run the
# command:   ./manage.py echo -p /data/echo -m /data/echo/echo-headers.yml

# you can try it out with the LLD compound library. Make sure to add plates
# with the barcode prefix BAF210901 and the correct dimensions to the
# experiment first. This can be done in the UI.

default_column_headers = {"source_plate_name": "Source Plate Name",
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
    "current_fluid_volume": "Current Fluid Volume", "percent_dms": "% DMSO"}


def convert_position_to_index(position: str, number_of_columns: int) -> int:
    match = re.match(r'([a-zA-Z]+)(\d+)', position)
    if match:
        letter = match.group(1)
        col = int(match.group(2))
        row = ord(letter.upper()) - 64
        index = row * number_of_columns + col
        return index


def read_csv_file(reader: csv.reader) -> List[Dict[str, str]]:
    """
    Reads a CSV file and returns a list of dictionaries representing each row in the file.
    """
    table_data = []
    column_headers = default_column_headers
    for row_index, table_row in enumerate(reader):
        if len(table_row) > 0 and table_row[0] == 'Source Plate Name':
            column_headers = table_row

        if row_index > 9:
            row_object = {}
            for cell_index, cell in enumerate(table_row):
                row_object[column_headers[cell_index]] = cell
            table_data.append(row_object)
    return table_data


def get_csv_files(path: str) -> list[dict[str, list[dict[str, str]] | str]]:
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
                                          "data": read_csv_file(csv_reader)})
    return files_data


# ./manage.py echo -p /data/echo -e exp
def parse_plate_data(mapping_data, file_path, headers):
    source_plate_barcode = mapping_data[0][headers['source_plate_barcode']]
    destination_plate_barcode = mapping_data[0][
        headers['destination_plate_barcode']]

    plates = Plate.objects.all()
    source_plate = plates.filter(barcode=source_plate_barcode).first()

    destination_plate = plates.filter(
        barcode=destination_plate_barcode).first()

    if not destination_plate:
        raise ValueError(
            f"Destination plate with barcode '{destination_plate_barcode}' "
            f"not found. Please add plates to the experiment.")
    number_of_columns = destination_plate.dimension.cols

    mapping_list = MappingList()
    for entry in mapping_data:
        if len(entry) > 0:
            log.debug(f"Mapping {entry[headers['source_well']]} to "
                      f"{entry[headers['destination_well']]} from"
                      f" {source_plate_barcode} to {destination_plate_barcode}")

            from_pos = convert_position_to_index(entry[headers['source_well']],
                                                 number_of_columns)
            to_pos = convert_position_to_index(
                entry[headers['destination_well']], number_of_columns)
            mapping = Mapping(from_pos=from_pos, to_pos=to_pos,
                              amount=float(entry[headers['actual_volume']]))
            mapping_list.add(mapping)
    mapping_success = source_plate.map(mapping_list, destination_plate)
    if mapping_success:
        file_content = open(file_path, 'rb').read()
        file_name = os.path.basename(file_path)
        plate_mapping = PlateMapping(source_plate=source_plate,
            target_plate=destination_plate)
        plate_mapping.mapping_file.save(file_name, ContentFile(file_content))
        log.info(f"Successfully mapped {source_plate_barcode} to "
                 f"{destination_plate_barcode}")



class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--path', '-p', type=str, required=True,
                            help='Path to the directory with the echo output')
        parser.add_argument('--experiment', '-e', type=str, required=False)
        parser.add_argument('--mapping-file', '-m', type=str,
            help='A yml file with the column headers, otherwise default '
                 'headers are used used')

    def handle(self, *args, **options):
        headers = default_column_headers
        headers_file = options.get('headers_file', None)
        if headers_file:
            try:
                with open(options.get('headers_file'), 'r') as file:
                    headers = yaml.safe_load(file)
            except FileNotFoundError:
                log.error(
                    f"The headers file '{headers_file}' could not be found.")
                return
            except yaml.YAMLError as e:
                log.error(f"Error parsing the YAML file '{headers_file}': {e}")
                return
        print('HEADERS: ', headers)

        try:
            path = options.get('path')
            data = get_csv_files(path)
            for mapping in data:
                parse_plate_data(mapping['data'], mapping['file_path'],
                                 headers)

        except Exception as ex:
            log.error(ex)
            traceback.print_exc()
