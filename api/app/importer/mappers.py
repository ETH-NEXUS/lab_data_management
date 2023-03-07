import csv
import json
import os
from glob import glob
from itertools import dropwhile
from tqdm import tqdm

import chardet
from core.mapping import Mapping, MappingList
from core.models import (
    BarcodeSpecification,
    Measurement,
    MeasurementFeature,
    MeasurementMetadata,
    Plate,
    PlateDimension,
    PlateMapping,
    Well,
)
from core.config import Config
from django.core.files import File
from friendlylog import colored_logger as log

from .helper import row_col_from_name


class BaseMapper:
    @staticmethod
    def get_files(_glob: str) -> list[str]:
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
                data = self.parse(filename, file, **kwargs)
            kwargs.update({"filename": filename})
            self.map(data, **kwargs)

    def parse(self, filename, file, **kwargs) -> dict:
        raise NotImplementedError

    def map(self, data: dict, **kwargs):
        raise NotImplementedError


class EchoMapper(BaseMapper):
    DEFAULT_COLUMNS = Config.current.importer.echo.default.columns

    def __fast_forward_to_header_row(self, file, headers):
        """
        Fast forward to header column determined
        by finding the first column name
        """
        pattern = list(headers.values())[0]
        return dropwhile(lambda line: pattern not in line, file)

    def parse(self, _, file, **kwargs):
        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        file = self.__fast_forward_to_header_row(file, headers)
        result = []
        reader = csv.DictReader(file, delimiter=",")
        reader = csv.DictReader(file, delimiter=",")

        for row in reader:
            # If there are None values in any of the following keys
            # or if there is a second header column we continue
            must_keys = (
                "source_plate_barcode",
                "source_plate_type",
                "source_well",
                "destination_plate_name",
                "destination_plate_barcode",
                "destination_well",
                "actual_volume",
            )
            if any(
                [
                    row[headers.get(key)] is None
                    or row[headers.get(key)] == headers.get(key)
                    for key in must_keys
                ]
            ):
                continue

            result.append(
                {new_key: row[old_key] for new_key, old_key in headers.items()}
            )
        return result

    def map(self, data: list[dict], **kwargs):
        def __debug(msg):
            if kwargs.get("debug"):
                log.debug(msg)

        # Plates cache by barcode
        plates = {}
        # MappingList cache by source plate barcode
        mapping_lists = {}
        # Process later queue
        queue = []

        with tqdm(
            desc="Processing mappings",
            unit="mappings",
            total=len(data),
        ) as mbar:
            for entry in data:
                source_plate_barcode = data[0]["source_plate_barcode"]
                destination_plate_name = data[0]["destination_plate_name"]
                destination_plate_barcode = data[0]["destination_plate_barcode"]

                if source_plate_barcode in plates:
                    source_plate = plates.get(source_plate_barcode)
                else:
                    try:
                        source_plate = Plate.objects.get(barcode=source_plate_barcode)
                        plates[source_plate_barcode] = source_plate
                    except Plate.DoesNotExist:
                        log.warning(
                            f"""Source plate with barcode {data[0]['source_plate_barcode']} does not exist.
                            I try again later..."""
                        )
                        queue.append(entry)
                        continue

                if destination_plate_barcode in plates:
                    destination_plate = plates.get(destination_plate_barcode)
                else:
                    try:
                        destination_plate = Plate.objects.get(
                            barcode=destination_plate_barcode
                        )
                    except Plate.DoesNotExist:
                        destination_plate = self.__create_plate_by_name_and_barcode(
                            destination_plate_name, destination_plate_barcode
                        )

                if source_plate_barcode in mapping_lists:
                    mapping_list = mapping_lists.get(source_plate_barcode)
                else:
                    mapping_list = MappingList(target=destination_plate)
                    mapping_lists[source_plate_barcode] = mapping_list

                source_well = entry["source_well"]
                destination_well = entry["destination_well"]
                __debug(
                    f"Mapping {source_plate.barcode}:{source_well} -> {destination_plate.barcode}:{destination_well}"
                )
                from_pos = source_plate.dimension.position(source_well)
                to_pos = destination_plate.dimension.position(destination_well)

                mapping = Mapping(
                    from_pos=from_pos,
                    to_pos=to_pos,
                    amount=float(entry["actual_volume"]),
                    status=entry["transfer_status"],
                )
                mapping_list.add(mapping)

                mbar.update(1)

        for source_plate_barcode, mapping_list in mapping_lists.items():
            source_plate = plates.get(source_plate_barcode)
            if source_plate.map(mapping_list, mapping_list.target):
                with open(kwargs["filename"], "rb") as file:
                    PlateMapping.objects.create(
                        source_plate=source_plate,
                        target_plate=mapping_list.target,
                        mapping_file=File(file, os.path.basename(file.name)),
                    )
                log.info(
                    f"Mapped {source_plate.barcode} -> {mapping_list.target.barcode}"
                )
            else:
                log.error(
                    f"Error mapping {source_plate.barcode} -> {mapping_list.target.barcode}"
                )

        # If there are entries queued because the source plate did not yet exist
        # we try map those again but only once.
        MAX_TRY_QUEUE = 3
        try_queue = kwargs.get("try_queue", 0)
        if try_queue < MAX_TRY_QUEUE and len(queue) > 0:
            kwargs.update({"try_queue": try_queue + 1})
            self.map(queue, **kwargs)

    def __create_plate_by_name_and_barcode(self, plate_name: str, barcode: str):
        try:
            barcode_prefix = barcode.split("_")[0]
            barcode_prefix = barcode.split("_")[0]
            barcode_specification = BarcodeSpecification.objects.get(
                prefix=barcode_prefix
            )
        except BarcodeSpecification.DoesNotExist:
            raise ValueError(f"No barcode specification found for {barcode}.")

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

    def map(self, data: dict, **kwargs):
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
