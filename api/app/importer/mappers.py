import csv
import json
import os
from glob import glob
from typing import List, Dict
from django.core.files.base import ContentFile
from django.core.exceptions import ObjectDoesNotExist
from friendlylog import colored_logger as log
import chardet

from core.mapping import Mapping, MappingList
from core.models import Plate, PlateMapping, Measurement, Well, \
    MeasurementMetadata, MeasurementFeature, MappingError, \
    BarcodeSpecification, PlateDimension


class BaseMapper:
    @staticmethod
    def get_files(_glob: str) -> List[str]:
        """Get all files in the given path that match the glob pattern"""
        return glob(_glob)

    def run(self, _glob, **kwargs):
        """Run the mapper"""
        for filename in self.get_files(_glob):
            log.info(f"Processing file {filename}...")
            with open(filename, 'rb') as file:
                encoding = chardet.detect(file.read()).get('encoding')
            with open(filename, 'r', encoding=encoding) as file:
                data = self.parse(filename, file, headers=kwargs['headers'])
            self.map(data, filename=filename)

    def parse(self, filename, file, **kwargs) -> Dict:
        raise NotImplementedError()

    def map(self, data: Dict, **kwargs):
        raise NotImplementedError()


class EchoMapper(BaseMapper):
    DEFAULT_COLUMNS = {"source_plate_barcode": "Source Plate Barcode",
        "source_plate_type": "Source Plate Type", "source_well": "Source Well",
        "destination_plate_name": "Destination Plate Name",
        "destination_plate_barcode": "Destination Plate Barcode",
        "destination_well": "Destination Well",
        "actual_volume": "Actual Volume"}

    def __fast_forward_to_header_row(self, file, headers):
        row = next(file)
        while 'DETAILS' not in row:  # tuple(headers.values())[0]
            row = next(file)
        return file

    def parse(self, filename, file, **kwargs):
        headers = kwargs.get('headers', EchoMapper.DEFAULT_COLUMNS)
        self.__fast_forward_to_header_row(file, headers)
        result = []
        reader = csv.DictReader(file, delimiter=',')

        for row in reader:
            result.append({new_key: row[old_key] for new_key, old_key in
                              headers.items()})
        file.close()
        return result

    def map(self, data: List[Dict], **kwargs):
        try:
            source_plate = Plate.objects.get(
                barcode=data[0]['source_plate_barcode'])
            if source_plate.dimension is None:
                source_plate.dimension = self.__find_plate_dimension(
                    data[0]['source_plate_name'])
                source_plate.save()
        except ObjectDoesNotExist:
            raise MappingError(f"Source plate with barcode does not "
                               f"exist.")
        try:
            destination_plate = Plate.objects.get(
                barcode=data[0]['destination_plate_barcode'])
        except ObjectDoesNotExist:
            destination_plate = self.__create_plate_by_name_and_barcode(
                plate_name=data[0]['destination_plate_name'],
                barcode=data[0]['destination_plate_barcode'])

        mapping_list = MappingList()

        for entry in data:
            print('ENTRY: ', entry)
            if len(entry) > 0:
                log.info(f"Mapping {entry['source_well']} to "
                         f"{entry['destination_well']} from"
                         f" {source_plate.barcode} to "
                         f"{destination_plate.barcode}")

                from_pos = source_plate.dimension.position(
                    entry['source_well'])
                to_pos = destination_plate.dimension.position(
                    entry['destination_well'])
                mapping = Mapping(from_pos=from_pos, to_pos=to_pos,
                                  amount=float(entry['actual_volume']))
                mapping_list.add(mapping)
        mapping_success = source_plate.map(mapping_list, destination_plate)
        if mapping_success:
            file_content = open(kwargs['filename'], 'rb').read()
            file_name = os.path.basename(kwargs['filename'])
            plate_mapping = PlateMapping(source_plate=source_plate,
                                         target_plate=destination_plate)
            plate_mapping.mapping_file.save(file_name,
                                            ContentFile(file_content))
            log.info(f"Successfully mapped {source_plate.barcode} to "
                     f"{destination_plate.barcode}")
        else:
            log.error(f"Failed to map {source_plate.barcode} to "
                      f"{destination_plate.barcode}")

    def __create_plate_by_name_and_barcode(self, plate_name: str, barcode:
    str):
        try:
            barcode_prefix = barcode.split('_')[0]
            barcode_specification = BarcodeSpecification.objects.get(
                prefix=barcode_prefix)
        except ObjectDoesNotExist:
            raise ValueError(f"No barcode specification found for {barcode}. "
                             f"Please add it in the user interface.")

        return Plate.objects.create(barcode=barcode,
                                    experiment=barcode_specification.experiment,
                                    dimension=self.__find_plate_dimension(
                                        plate_name))


    def __find_plate_dimension(self, plate_name: str):
        # default and the most common dimension
        plate_dimension_name = 'dim_384_16x24'
        try:
            if '96' in plate_name:
                plate_dimension_name = 'dim_96_8x12'
            elif '1536' in plate_name:
                plate_dimension_name = 'dim_1536_32x48'
            plate_dimension = PlateDimension.objects.get(
                name=plate_dimension_name)
            return plate_dimension
        except ObjectDoesNotExist:
            raise ValueError(
                f"No plate dimension found for name {plate_dimension_name}. "
                f"Please add it in the user interface.")


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
                                                                     index][0],
                                                                 feature=
                                                                 metadata_list[
                                                                     index][1])
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
                    measurement_feature_name = metadata[
                        index + 1].strip().lstrip()
        measurement_feature = MeasurementFeature.objects.get_or_create(
            name=measurement_feature_name)[0]

        for item in dict_list:
            measurement_metadata = MeasurementMetadata.objects.create(
                data=json.dumps(item))
            log.debug(f"Successfully created metadata")
            metadata_objects_list.append(
                (measurement_metadata, measurement_feature))
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
