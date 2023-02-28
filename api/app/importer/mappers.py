import csv
import json
import os
from glob import glob
from typing import List, Dict
from django.core.files.base import ContentFile
from friendlylog import colored_logger as log

from core.mapping import Mapping, MappingList
from core.models import Plate, PlateMapping, Measurement, Well, \
    MeasurementMetadata, MeasurementFeature

class BaseMapper:
    @staticmethod
    def get_files(self, _glob:str) -> List[str]:
        """Get all files in the given path that match the glob pattern"""
        return glob(_glob)

class EchoMapper(BaseMapper):
    DEFAULT_COLUMN_HEADERS = {
        "source_plate_name": "Source Plate Name",
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
        "percent_dms": "% DMSO"
    }

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
                                              "data":
                                                  EchoMapper.read_csv_echo_file(
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


class MeasurementMapper:

    @staticmethod
    def get_measurement_files(path: str):
        files_data = []
        for root, dirs, files in os.walk(path, topdown=False):
            for file in files:
                if file.endswith(".asc"):
                    file_path = os.path.join(root, file)
                    barcode_delimiter_index = file.find('_')
                    barcode = file[barcode_delimiter_index + 1: -4]
                    files_data.append(
                        {"barcode": barcode, "file_path": file_path,
                            "measurement_data": MeasurementMapper.read_measurement_file(
                                file_path)})
                    log.info(f"Read file  {file_path} ")

        return files_data

    @staticmethod
    def read_measurement_file(file_path: str) -> Dict[str, List[str]]:
        with open(file_path, 'r', encoding="iso-8859-1") as measurement_file:
            measurement_data = measurement_file.readlines()
            metadata_start_index = MeasurementMapper.find_metadata_start_index(
                measurement_data)
            metadata_list = MeasurementMapper.parse_measurement_metadata(
                measurement_data[metadata_start_index:])
            return {"values": measurement_data[: metadata_start_index],
                "metadata": metadata_list}

    @staticmethod
    def parse_measurement_data(measurement_data: Dict[str, List[str]],
                               barcode: str) -> None:

        plate = Plate.objects.get(barcode=barcode)

        if not plate:
            raise ValueError(f"Plate with barcode '{barcode}' "
                             f"not found. Please add plates to the experiment.")
        number_of_columns = plate.dimension.cols

        # We drop the first line if it contains 'A1' because it is the header
        if 'A1' in measurement_data['values'][0].strip().split('\t'):
            first_line = measurement_data['values'][0].strip().split('\t')
        else:
            first_line = measurement_data['values'][1].strip().split('\t')
        indices = MeasurementMapper.find_indices(first_line)

        metadata_list = measurement_data['metadata']

        for line in measurement_data['values']:
            try:
                line_list = line.strip().split('\t')

                well_position_str = line_list[indices[0]]
                well_position = Plate.convert_position_to_index(
                    well_position_str.strip().lstrip(), number_of_columns)
                identifier = line_list[indices[1]]
                values = line_list[2:]
                well = Well.objects.filter(position=well_position,
                                           plate=plate).first()

                if well:
                    for index, value in enumerate(values):
                        pass
                        measurement = Measurement.objects.create(well=well,
                                                                 value=value,
                                                                 identifier=identifier,
                                                                 meta=
                                                                 metadata_list[
                                                                     index][
                                                                     0],
                                                                 feature=metadata_list[index][1])
                        measurement.save()
                        log.debug(f"Succesfully created measurement for well "
                                  f"{well_position} with value {value}")
            except (ValueError, Well.DoesNotExist) as e:
                log.error(f"Error processing line: {line.strip()}. {str(e)}")

    @staticmethod
    def parse_measurement_metadata(metadata):
        """
        Creates metadata objects for every value of measurement.
        """
        metadata_objects_list = []
        dict_list = [{}]
        measurement_feature_name = 'unknown'
        for index, line in enumerate(metadata):
            if ":" in line:
                key, value = line.split(':', 1)
                key = key.strip().lstrip()
                value = value.strip().lstrip()
                if key in dict_list[-1]:
                    dict_list.append({})
                dict_list[-1][key] = value
                if key == 'Range':
                    measurement_feature_name = metadata[index + 1].strip().lstrip()
        measurement_feature = MeasurementFeature.objects.get_or_create(name=measurement_feature_name)[0]

        for item in dict_list:
            measurement_metadata = MeasurementMetadata.objects.create(data=json.dumps(item))
            log.debug(f"Successfully created metadata")
            metadata_objects_list.append((measurement_metadata, measurement_feature))
        return metadata_objects_list

    @staticmethod
    def find_indices(line_list: list[str]) -> tuple[int, int]:
        """
        find, which column has the identifier and which has the well (we
        do it by checking the first line and finding the index of 'A1'
        we can not use regular expressions because there are identifiers
        like 'NC1' which is the same pattern as by well names
        """
        well_index = 0
        identifier_index = 1
        for index, item in enumerate(line_list):
            # If the item is a numeric character or string, skip it, because
            # it is a value
            if item.isnumeric():
                continue
            if item == 'A1':
                well_index = index
            else:
                identifier_index = index
        return well_index, identifier_index

    @staticmethod
    def find_metadata_start_index(measurement_data: list[str]) -> int:
        """
        in the lines of the measurement file, finds the index of the line where
        metadata starts
        """
        index = -1
        for index, line in enumerate(measurement_data):
            if line.startswith('Date of measurement'):
                index = index
                break
        return index