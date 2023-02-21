from django.core.management.base import BaseCommand
from friendlylog import colored_logger as log
import traceback
import csv
import os
from typing import List, Dict
from core.models import Experiment, Plate, PlateMapping
from core.mapping import Mapping, MappingList, PositionMapper
import re

default_column_headers = ['Source Plate Name', 'Source Plate Barcode',
                          'Source Plate Type', 'Source Well',
                          'Source Concentration', 'Source Concentration Units',
                          'Destination Plate Name',
                          'Destination Plate Barcode', 'Destination Well',
                          'Destination Concentration',
                          'Destination Concentration Units', 'Compound Name',
                          'Transfer Volume', 'Actual Volume',
                          'Transfer Status', 'Current Fluid Height',
                          'Current Fluid Volume', '% DMSO']

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


def get_csv_files(path: str) -> List[List[Dict[str, str]]]:
    """
    Recursively searches the given directory for CSV files with "transfer" in their name,
    reads each file using the read_csv_file function, and returns a list of tables
    representing the data in the CSV files.
    """
    files_data = []
    for root, dirs, files in os.walk(path, topdown=False):
        for file in files:
            if file.endswith(".csv") and 'transfer' in file:
                with open(os.path.join(root, file), 'r') as csv_file:
                    log.debug(f"Reading file {file}")
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    files_data.append(read_csv_file(csv_reader))
    return files_data


# ./manage.py echo -p /data/echo -e exp
def parse_plate_data(mapping_data):
    source_plate_barcode = mapping_data[0]['Source Plate Barcode']
    destination_plate_barcode = mapping_data[0]['Destination Plate Barcode']

    plates = Plate.objects.all()
    source_plate = plates.filter(barcode=source_plate_barcode).first()
    destination_plate = plates.filter(
        barcode=destination_plate_barcode).first()
    number_of_columns = destination_plate.dimension.cols


    mapping_list = MappingList()
    for entry in mapping_data:
        if len(entry) > 0:
            log.debug(
                f"Mapping {entry['Source Well']} to {entry['Destination Well']} from {source_plate_barcode} to {destination_plate_barcode}")


            from_pos = convert_position_to_index(entry['Source Well'], number_of_columns)
            to_pos = convert_position_to_index(entry['Destination Well'], number_of_columns)

            mapping = Mapping(from_pos=from_pos,
                              to_pos=to_pos,
                              amount=float(entry['Actual Volume']))
            mapping_list.add(mapping)
    mapping_success = source_plate.map(mapping_list, destination_plate)
    if mapping_success:
        PlateMapping.objects.create(source_plate=source_plate,
                                    target_plate=destination_plate)
        log.info(f"Successfully mapped {source_plate_barcode} to "
                 f"{destination_plate_barcode}")



class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--path', '-p', type=str, required=True,
                            help='Path to the directory with the echo output')
        parser.add_argument('--experiment', '-e', type=str, required=False)

    def handle(self, *args, **options):
        try:
            path = options.get('path')
            data = get_csv_files(path)
            for mapping in data:
                parse_plate_data(mapping)

        except Exception as ex:
            log.error(ex)
            traceback.print_exc()
