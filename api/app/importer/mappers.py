import csv
import json
import os
from glob import glob
from typing import Dict, List

import chardet
from core.mapping import Mapping, MappingList
from core.models import (
    BarcodeSpecification,
    MappingError,
    Measurement,
    MeasurementFeature,
    MeasurementMetadata,
    Plate,
    PlateDimension,
    PlateMapping,
    Well,
)
from django.core.exceptions import ObjectDoesNotExist
from django.core.files.base import ContentFile
from friendlylog import colored_logger as log

from .helper import row_col_from_name


class BaseMapper:
    @staticmethod
    def get_files(_glob: str) -> List[str]:
        """Get all files in the given path that match the glob pattern"""
        return glob(_glob)

    def run(self, _glob, **kwargs):
        """Run the mapper"""
        for filename in self.get_files(_glob):
            log.info(f"Processing file {filename}...")
            with open(filename, "rb") as file:
                encoding = chardet.detect(file.read()).get("encoding")
            with open(filename, "r", encoding=encoding) as file:
            with open(filename, "rb") as file:
                encoding = chardet.detect(file.read()).get("encoding")
            with open(filename, "r", encoding=encoding) as file:
                data = self.parse(
                    filename,
                    file,
                    headers=kwargs["headers"] if "headers" in kwargs else None,
                    filename,
                    file,
                    headers=kwargs["headers"] if "headers" in kwargs else None,
                )
            self.map(data, filename=filename)

    def parse(self, filename, file, **kwargs) -> Dict:
        raise NotImplementedError

    def map(self, data: Dict, **kwargs):
        raise NotImplementedError


class EchoMapper(BaseMapper):
    DEFAULT_COLUMNS = {
        "source_plate_barcode": "Source Plate Barcode",
        "source_plate_type": "Source Plate Type",
        "source_well": "Source Well",
        "source_plate_type": "Source Plate Type",
        "source_well": "Source Well",
        "destination_plate_name": "Destination Plate Name",
        "destination_plate_barcode": "Destination Plate Barcode",
        "destination_well": "Destination Well",
        "actual_volume": "Actual Volume",
        "actual_volume": "Actual Volume",
    }

    def __fast_forward_to_header_row(self, file, headers):
        row = next(file)
        while "DETAILS" not in row:  # tuple(headers.values())[0]
        while "DETAILS" not in row:  # tuple(headers.values())[0]
            row = next(file)
        return file

    def parse(self, filename, file, **kwargs):
        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        self.__fast_forward_to_header_row(file, headers)
        result = []
        reader = csv.DictReader(file, delimiter=",")
        reader = csv.DictReader(file, delimiter=",")

        for row in reader:
            result.append(
                {new_key: row[old_key] for new_key, old_key in headers.items()}
            )
        file.close()
        return result

    def map(self, data: List[Dict], **kwargs):
        try:
            source_plate = Plate.objects.get(barcode=data[0]["source_plate_barcode"])
            source_plate = Plate.objects.get(barcode=data[0]["source_plate_barcode"])
            if source_plate.dimension is None:
                source_plate.dimension = self.__get_plate_dimension(
                    data[0]["source_plate_name"]
                )
                source_plate.save()
        except Plate.DoesNotExist:
            raise MappingError(
                f"Source plate with barcode {data[0]['source_plate_barcode']} does not exist."
            )
        try:
            destination_plate = Plate.objects.get(
                barcode=data[0]["destination_plate_barcode"]
                barcode=data[0]["destination_plate_barcode"]
            )
        except ObjectDoesNotExist:
            destination_plate = self.__create_plate_by_name_and_barcode(
                plate_name=data[0]["destination_plate_name"],
                barcode=data[0]["destination_plate_barcode"],
                plate_name=data[0]["destination_plate_name"],
                barcode=data[0]["destination_plate_barcode"],
            )

        mapping_list = MappingList()

        for entry in data:
            if len(entry) > 0:
                log.info(
                    f"Mapping {entry['source_well']} to "
                    f"{entry['destination_well']} from"
                    f" {source_plate.barcode} to "
                    f"{destination_plate.barcode}"
                )
                from_pos = source_plate.dimension.position(entry["source_well"])
                to_pos = destination_plate.dimension.position(entry["destination_well"])
                mapping = Mapping(
                    from_pos=from_pos,
                    to_pos=to_pos,
                    amount=float(entry["actual_volume"]),
                    from_pos=from_pos,
                    to_pos=to_pos,
                    amount=float(entry["actual_volume"]),
                )
                mapping_list.add(mapping)
        mapping_success = source_plate.map(mapping_list, destination_plate)
        if mapping_success:
            file_content = open(kwargs["filename"], "rb").read()
            file_name = os.path.basename(kwargs["filename"])
            file_content = open(kwargs["filename"], "rb").read()
            file_name = os.path.basename(kwargs["filename"])
            plate_mapping = PlateMapping(
                source_plate=source_plate, target_plate=destination_plate
            )
            plate_mapping.mapping_file.save(file_name, ContentFile(file_content))
            plate_mapping.mapping_file.save(file_name, ContentFile(file_content))
            log.info(
                f"Successfully mapped {source_plate.barcode} to "
                f"{destination_plate.barcode}"
            )
        else:
            log.error(
                f"Failed to map {source_plate.barcode} to "
                f"{destination_plate.barcode}"
            )

    def __create_plate_by_name_and_barcode(self, plate_name: str, barcode: str):
    def __create_plate_by_name_and_barcode(self, plate_name: str, barcode: str):
        try:
            barcode_prefix = barcode.split("_")[0]
            barcode_prefix = barcode.split("_")[0]
            barcode_specification = BarcodeSpecification.objects.get(
                prefix=barcode_prefix
            )
        except ObjectDoesNotExist:
            raise ValueError(
                f"No barcode specification found for {barcode}. "
                f"Please add it in the user interface."
            )

        return Plate.objects.create(
            barcode=barcode,
            experiment=barcode_specification.experiment,
            dimension=self.__get_plate_dimension(plate_name),
        )

    def __get_plate_dimension(self, plate_name: str):
        try:
            rows, cols = row_col_from_name(plate_name)
            plate_dimension = PlateDimension.objects.get(rows=rows, cols=cols)
            return plate_dimension
        except ValueError:
            raise
        except PlateDimension.DoesNotExist:
            raise ValueError(f"No plate dimension found: {rows}x{cols}")


class MeasurementMapper(BaseMapper):
    def parse(self, filename, file, **kwargs):
        log.info(f"Read file  {filename} ")
        barcode_delimiter_index = filename.find("_")
        barcode = filename[barcode_delimiter_index + 1 : -4]
        barcode_delimiter_index = filename.find("_")
        barcode = filename[barcode_delimiter_index + 1 : -4]
        data = file.readlines()
        metadata_start_index = self.__find_metadata_start_index(data)
        metadata_list = self.__parse_measurement_metadata(data[metadata_start_index:])
        metadata_list = self.__parse_measurement_metadata(data[metadata_start_index:])

        return {
            "barcode": barcode,
            "filename": filename,
            "measurement_data": {
            "barcode": barcode,
            "filename": filename,
            "measurement_data": {
                "values": data[:metadata_start_index],
                "metadata": metadata_list,
            },
                "metadata": metadata_list,
            },
        }

    def map(self, data: Dict, **kwargs):
        try:
            plate = Plate.objects.get(barcode=data["barcode"])
        except Plate.DoesNotExist:
            raise ValueError(
                f"Plate with barcode {data['barcode']} not found. "
                f"Please add plates to the experiment."
            )
        else:
            values = data["measurement_data"]["values"]

            # We drop the first line if it contains 'A1' because it is the
            # header
            if "A1" in values[0].strip().split("\t"):
                first_line = values[0].strip().split("\t")
            else:
                first_line = values[1].strip().split("\t")
                first_line = values[1].strip().split("\t")
            indices = self.__find_indices(first_line)
            metadata_list = data["measurement_data"]["metadata"]
            for line in values:
                try:
                    line_list = line.strip().split("\t")
                    line_list = line.strip().split("\t")

                    well_position_str = line_list[indices[0]]
                    well_position = plate.dimension.position(
                        well_position_str.strip().lstrip()
                    )
                        well_position_str.strip().lstrip()
                    )
                    identifier = line_list[indices[1]]
                    values = line_list[2:]
                    well = Well.objects.filter(
                        position=well_position, plate=plate
                    ).first()

                    if well:
                        for index, value in enumerate(values):
                            measurement = Measurement.objects.create(
                                well=well,
                                value=value,
                                identifier=identifier,
                                well=well,
                                value=value,
                                identifier=identifier,
                                meta=metadata_list[index][0],
                                feature=metadata_list[index][1],
                                feature=metadata_list[index][1],
                            )
                            measurement.save()
                            log.debug(
                                f"Succesfully created measurement for well "
                                f"{well_position} with value {value}"
                            )
                except (ValueError, Well.DoesNotExist) as e:
                    log.error(f"Error processing line: {line.strip()}. {str(e)}")
                    log.error(f"Error processing line: {line.strip()}. {str(e)}")

    def __parse_measurement_metadata(self, metadata):
        """
        Creates metadata objects for every value of measurement.
        """
        metadata_objects_list = []
        dict_list = [{}]
        measurement_feature_name = "unknown"
        measurement_feature_name = "unknown"
        for index, line in enumerate(metadata):
            if ":" in line:
                key, value = line.split(":", 1)
                key, value = line.split(":", 1)
                key = key.strip().lstrip()
                value = value.strip().lstrip()
                if key in dict_list[-1]:
                    dict_list.append({})
                dict_list[-1][key] = value
                if key == "Range":
                    measurement_feature_name = metadata[index + 1].strip().lstrip()
                if key == "Range":
                    measurement_feature_name = metadata[index + 1].strip().lstrip()
        measurement_feature = MeasurementFeature.objects.get_or_create(
            name=measurement_feature_name
        )[0]

        for item in dict_list:
            measurement_metadata = MeasurementMetadata.objects.create(
                data=json.dumps(item)
            )
            log.debug("Successfully created metadata")
            metadata_objects_list.append((measurement_metadata, measurement_feature))
        return metadata_objects_list

    def __find_indices(self, line_list: list[str]) -> tuple[int, int]:
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
            if item == "A1":
            if item == "A1":
                well_index = index
            else:
                identifier_index = index
        return well_index, identifier_index

    def __find_metadata_start_index(self, measurement_data: list[str]) -> int:
        """
        in the lines of the measurement file, finds the index of the line where
        metadata starts
        """
        index = -1
        for index, line in enumerate(measurement_data):
            if line.startswith("Date of measurement"):
            if line.startswith("Date of measurement"):
                index = index
                break
        return index