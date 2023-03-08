import csv
import os
import re
from io import TextIOWrapper
from glob import glob
from itertools import dropwhile
from tqdm import tqdm
from datetime import datetime as dt
from contextlib import redirect_stderr

from chardet.universaldetector import UniversalDetector
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
    MappingError,
)
from core.config import Config
from django.core.files import File
from django.utils import timezone
from friendlylog import colored_logger as log

from .helper import row_col_from_name

import logging

logging.getLogger("chardet.charsetprober").setLevel(logging.INFO)


class BaseMapper:
    @staticmethod
    def get_files(_glob: str) -> list[str]:
        """Get all files in the given path that match the glob pattern"""
        return glob(_glob)

    def run(self, _glob, **kwargs):
        """Run the mapper"""
        for filename in self.get_files(_glob):
            log.info(f"Processing file {filename}...")
            with redirect_stderr(None):
                detector = UniversalDetector()
                with open(filename, "rb") as file:
                    for line in file:
                        detector.feed(line)
                        if detector.done:
                            break
                    detector.close()
                encoding = detector.result.get("encoding")
            with open(filename, "r", encoding=encoding) as file:
                ret = self.parse(file, **kwargs)
                if isinstance(ret, tuple):
                    data = ret[0]
                    kwargs.update(ret[1])
                else:
                    data = ret
            kwargs.update({"filename": filename})
            self.map(data, **kwargs)

    def parse(self, file: TextIOWrapper, **kwargs) -> list[dict]:
        raise NotImplementedError

    def map(self, data: list[dict], **kwargs) -> None:
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

    def parse(self, file: TextIOWrapper, **kwargs) -> list[dict]:
        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        file = self.__fast_forward_to_header_row(file, headers)
        results = []
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

            results.append(
                {new_key: row[old_key] for new_key, old_key in headers.items()}
            )
        return results

    def map(self, data: list[dict], **kwargs) -> None:
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
                source_plate_barcode = entry["source_plate_barcode"]
                destination_plate_name = entry["destination_plate_name"]
                destination_plate_barcode = entry["destination_plate_barcode"]

                if source_plate_barcode in plates:
                    source_plate = plates.get(source_plate_barcode)
                else:
                    try:
                        source_plate = Plate.objects.get(barcode=source_plate_barcode)
                        plates[source_plate_barcode] = source_plate
                    except Plate.DoesNotExist:
                        log.warning(
                            f"""Source plate with barcode {entry['source_plate_barcode']} does not exist.
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


class M1000Mapper(BaseMapper):
    RE_FILENAME = r"(?P<date>[0-9]+)-(?P<time>[0-9]+)_(?P<barcode>[^\.]+)\.asc"
    RE_POS = r"^(?P<pos>[A-Z]+[0-9]+)\t(?P<rest>.+)\t$"
    RE_NUM = r"[0-9\.]+"
    RE_DATE_OF_MEASUREMENT = r"^Date of measurement: (?P<date>[^\/]+)\/Time of measurement: (?P<time>[0-9:]+)$"
    RE_PLATE_DESCRIPTION = r"^Plate Description: (?P<description>.+)$"
    RE_META_DATA = r"^    (?P<key>[^:]+): (?P<value>.+)$"

    def parse(self, file: TextIOWrapper, **kwargs) -> list[dict]:
        match = re.match(self.RE_FILENAME, os.path.basename(file.name))
        if match:
            barcode = match.group("barcode")
        else:
            raise MappingError(
                f"File name {os.path.basename(file.name)} does not much conventions."
            )

        results = []
        measurement_date = timezone.now()
        plate_description = None
        meta_data = [{}]
        meta_data_idx = 0
        for line in file:
            if match := re.match(self.RE_POS, line):
                parts = match.group("rest").split("\t")
                values = []
                identifier = None
                for part in parts:
                    if re.match(self.RE_NUM, part):
                        values.append(float(part))
                    else:
                        identifier = part
                results.append(
                    {
                        "position": match.group("pos"),
                        "values": values,
                        "identifier": identifier,
                    }
                )
            elif match := re.match(self.RE_DATE_OF_MEASUREMENT, line):
                measurement_date = dt.strptime(
                    f"{match.group('date')} {match.group('time')}", "%Y-%m-%d %H:%M:%S"
                )
            elif match := re.match(self.RE_PLATE_DESCRIPTION, line):
                plate_description = match.group("description")
            elif match := re.match(self.RE_META_DATA, line):
                if match.group("key") in meta_data[meta_data_idx]:
                    meta_data_idx += 1
                    meta_data.append({})
                meta_data[meta_data_idx][match.group("key")] = match.group("value")
            else:
                pass

        kwargs.update(
            {
                "barcode": barcode,
                "measurement_date": measurement_date,
                "plate_description": plate_description,
                "meta_data": meta_data,
            }
        )
        return results, kwargs

    def map(self, data: list[dict], **kwargs) -> None:
        barcode = kwargs.get("barcode")
        try:
            plate = Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:
            raise ValueError(f"Plate with barcode {barcode} does not exist.")

        with tqdm(
            desc="Processing measurements",
            unit="measurement",
            total=len(data),
        ) as mbar:
            for entry in data:
                position = plate.dimension.position(entry.get("position"))
                well = plate.well_at(position)
                if not well:
                    well = Well.objects.create(plate=plate, position=position)

                for idx, value in enumerate(entry.get("values")):
                    feature, _ = MeasurementFeature.objects.get_or_create(
                        abbrev=kwargs.get("meta_data")[idx].get("Label")
                    )
                    metadata, _ = MeasurementMetadata.objects.get_or_create(
                        data=kwargs.get("meta_data")[idx]
                    )
                    # log.debug(f"Add measurement {plate.barcode}{well.hr_position}: {value}")
                    Measurement.objects.update_or_create(
                        well=well,
                        feature=feature,
                        defaults={
                            "value": value,
                            "identifier": entry.get("identifier"),
                            "meta": metadata,
                        },
                    )
                mbar.update(1)

            # try:
            #     line_list = line.strip().split("\t")

            #     well_position_str = line_list[indices[0]]
            #     well_position = plate.dimension.position(
            #         well_position_str.strip().lstrip()
            #     )
            #     identifier = line_list[indices[1]]
            #     values = line_list[2:]
            #     well = Well.objects.filter(position=well_position, plate=plate).first()

            #     if well:
            #         for index, value in enumerate(values):
            #             measurement = Measurement.objects.create(
            #                 well=well,
            #                 value=value,
            #                 identifier=identifier,
            #                 meta=metadata_list[index][0],
            #                 feature=metadata_list[index][1],
            #             )
            #             measurement.save()
            #             log.debug(
            #                 f"Succesfully created measurement for well "
            #                 f"{well_position} with value {value}"
            #             )
            # except (ValueError, Well.DoesNotExist) as e:
            #     log.error(f"Error processing line: {line.strip()}. {str(e)}")

    # def __parse_measurement_metadata(self, metadata):
    #     """
    #     Creates metadata objects for every value of measurement.
    #     """
    #     metadata_objects_list = []
    #     dict_list = [{}]
    #     measurement_feature_name = "unknown"
    #     for index, line in enumerate(metadata):
    #         if ":" in line:
    #             key, value = line.split(":", 1)
    #             key = key.strip().lstrip()
    #             value = value.strip().lstrip()
    #             if key in dict_list[-1]:
    #                 dict_list.append({})
    #             dict_list[-1][key] = value
    #             if key == "Range":
    #                 measurement_feature_name = metadata[index + 1].strip().lstrip()
    #     measurement_feature = MeasurementFeature.objects.get_or_create(
    #         name=measurement_feature_name
    #     )[0]

    #     for item in dict_list:
    #         measurement_metadata = MeasurementMetadata.objects.create(
    #             data=json.dumps(item)
    #         )
    #         log.debug("Successfully created metadata")
    #         metadata_objects_list.append((measurement_metadata, measurement_feature))
    #     return metadata_objects_list

    # def __find_indices(self, line_list: list[str]) -> tuple[int, int]:
    #     """
    #     find, which column has the identifier and which has the well (we
    #     do it by checking the first line and finding the index of 'A1'
    #     we can not use regular expressions because there are identifiers
    #     like 'NC1' which is the same pattern as by well names
    #     """
    #     well_index = 0
    #     identifier_index = 1
    #     for index, item in enumerate(line_list):
    #         # If the item is a numeric character or string, skip it, because
    #         # it is a value
    #         if item.isnumeric():
    #             continue
    #         if item == "A1":
    #             well_index = index
    #         else:
    #             identifier_index = index
    #     return well_index, identifier_index

    # def __find_metadata_start_index(self, measurement_data: list[str]) -> int:
    #     """
    #     in the lines of the measurement file, finds the index of the line where
    #     metadata starts
    #     """
    #     index = -1
    #     for index, line in enumerate(measurement_data):
    #         if line.startswith("Date of measurement"):
    #             index = index
    #             break
    #     return index
