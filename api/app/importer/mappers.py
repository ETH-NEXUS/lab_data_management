import csv
import os
import re
import pandas as pd
from io import TextIOWrapper
from glob import glob
from itertools import dropwhile
from tqdm import tqdm
from datetime import datetime as dt
from contextlib import redirect_stderr
from importer.helper import message

from chardet.universaldetector import UniversalDetector
from core.mapping import Mapping, MappingList
from core.models import (
    BarcodeSpecification,
    Measurement,
    MeasurementFeature,
    Plate,
    PlateDimension,
    PlateMapping,
    Well,
    MappingError,
    MeasurementAssignment,
    Experiment,
    PlateDetail,
    WellDetail,
    ExperimentDetail,
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
        return glob(_glob, recursive=True)

    def run(self, _glob, **kwargs):
        """Run the mapper"""
        for filename in self.get_files(_glob):
            message(
                f"Processing file {filename}...", "info", kwargs.get("room_name", None)
            )
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

        message(
            "Refreshing materialized views...", "info", kwargs.get("room_name", None)
        )
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

    def parse(self, file: TextIOWrapper, **kwargs):
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
                message(msg, "debug", kwargs.get("room_name", None))

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
                destination_plate_type = entry["destination_plate_type"]

                destination_plate_barcode = entry["destination_plate_barcode"]

                if source_plate_barcode in plates:
                    source_plate = plates.get(source_plate_barcode)
                else:
                    try:
                        source_plate = Plate.objects.get(barcode=source_plate_barcode)
                        plates[source_plate_barcode] = source_plate
                    except Plate.DoesNotExist:
                        message(
                            f"""Source plate with barcode {entry['source_plate_barcode']} does not exist.
                            I try again later...""",
                            "warning",
                            kwargs.get("room_name", None),
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
                        message(
                            f"Creating destination plate {destination_plate_name}, {destination_plate_type}"
                        )
                        destination_plate = self.__create_plate_by_name_and_barcode(
                            destination_plate_name,
                            destination_plate_type,
                            destination_plate_barcode,
                            **kwargs,
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
                message(
                    f"Mapped {source_plate.barcode} -> {mapping_list.target.barcode}",
                    "info",
                    kwargs.get("room_name", None),
                )
            else:
                message(
                    f"Error mapping {source_plate.barcode} -> {mapping_list.target.barcode}",
                    "error",
                    kwargs.get("room_name", None),
                )

        # If there are entries queued because the source plate did not yet exist
        # we try map those again but only once.
        MAX_TRY_QUEUE = 3
        try_queue = kwargs.get("try_queue", 0)
        if try_queue < MAX_TRY_QUEUE and len(queue) > 0:
            kwargs.update({"try_queue": try_queue + 1})
            self.map(queue, **kwargs)

    def __create_plate_by_name_and_barcode(
        self, plate_name: str, plate_type: str, barcode: str, **kwargs
    ):
        try:
            barcode_prefix = barcode.split("_")[0]
            barcode_specification = BarcodeSpecification.objects.get(
                prefix=barcode_prefix
            )
        except BarcodeSpecification.DoesNotExist:
            if kwargs.get("experiment_name"):
                message(
                    f"No barcode specification found for {barcode}. Creating it.",
                    "warning",
                    kwargs.get("room_name", None),
                )

                barcode_specification, _ = BarcodeSpecification.objects.get_or_create(
                    prefix=barcode.split("_")[0],
                    sides=["North"],
                    number_of_plates=4,
                    experiment=Experiment.objects.get(
                        name=kwargs.get("experiment_name")
                    ),
                )
            else:
                message(
                    f"No barcode specification found for {barcode} and no experiment name is provided. Please provide the experiment name in order to create the missing barcode specifications.",
                    "error",
                    kwargs.get("room_name", None),
                )
                raise ValueError(
                    f"No barcode specification found for {barcode} and no experiment name is "
                    f"provided."
                    f" Please provide the experiment name in order to create the missing barcode "
                    f"specifications."
                )

        return Plate.objects.create(
            barcode=barcode,
            experiment=barcode_specification.experiment,
            dimension=self.__get_plate_dimension(
                plate_name, plate_type, kwargs.get("room_name")
            ),
        )

    def __get_plate_dimension(self, plate_name: str, plate_type: str, room_name):
        try:
            rows, cols = row_col_from_name(f"{plate_type} {plate_name}")
            plate_dimension = PlateDimension.objects.get(rows=rows, cols=cols)
            return plate_dimension
        except ValueError:
            raise
        except PlateDimension.DoesNotExist:
            message(f"No plate dimension found: {rows}x{cols}", "error", room_name)
            raise ValueError(f"No plate dimension found: {rows}x{cols}")


class M1000Mapper(BaseMapper):
    RE_FILENAME = r"(?P<date>[0-9]+)-(?P<time>[0-9]+)_(?P<barcode>[^\.]+)\.asc"

    RE_POS = r"^[A-Z]+[0-9]+$"
    RE_ID = r"^[^_]+_[^_]+$"
    RE_NUM = r"^[0-9\.]+$"
    RE_TAB = r"[\t\s]+"

    RE_DATE_OF_MEASUREMENT = r"^Date of measurement: (?P<date>[^\/]+)\/Time of measurement: (?P<time>[0-9:]+)$"
    RE_PLATE_DESCRIPTION = r"^Plate Description: (?P<description>.+)$"
    RE_META_DATA = r"^    (?P<key>[^:]+): (?P<value>.+)$"

    def determine_indexes(self, file: TextIOWrapper):
        """
        Determines the indexes for position and identifier.
        The rest of the columns are values.

        The strategy is to go through the lines of the file
        and stop as soon as we find an unambiguous line.
        """
        for line in file:
            parts = re.split(self.RE_TAB, line.strip())
            unambiguity = []
            pos_index = None
            id_index = None
            for idx, part in enumerate(parts):
                if match_pos := re.match(self.RE_POS, part):
                    pos_index = idx
                if match_id := re.match(self.RE_ID, part):
                    id_index = idx
                # if this is an unambiguous value
                unambiguity.append([match_pos, match_id])
            # if the line is an unambiguous line
            df = pd.DataFrame(unambiguity)
            if df.count().apply(lambda x: x == 1).all():
                file.seek(0)
                return pos_index, id_index
        file.seek(0)
        raise MappingError(f"File has not the desired format: {file.name}")

    def parse(self, file: TextIOWrapper, **kwargs) -> list[dict]:
        def __debug(msg):
            if kwargs.get("debug"):
                message(msg, "debug", kwargs.get("room_name", None))

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
        pos_index, id_index = self.determine_indexes(file)
        for line in file:
            parts = re.split(self.RE_TAB, line)
            # For a value line we expect at least 3 parts (pos, id, value)
            # and the position part should match the position regex
            if len(parts) >= 3 and re.match(self.RE_POS, parts[pos_index]):
                position = parts[pos_index]
                identifier = parts[id_index]
                result = {
                    "position": position,
                    "identifier": identifier,
                    "values": [
                        part
                        for idx, part in enumerate(parts)
                        if idx not in [pos_index, id_index]
                        and re.match(self.RE_NUM, part)
                    ],
                }
                # If there are no values in this line we ignore it
                if len(result["values"]) > 0:
                    __debug(f"result: {result}")
                    results.append(result)
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
                "meta_data": list(reversed(meta_data)),  # the order of the
                # metadata in the file is reversed compared to the order of the values
            }
        )
        return results, kwargs

    def apply_evaluation_formula(self, formula, plate, well, entry, **kwargs):
        result = 0
        for idx, value in enumerate(entry.get("values")):
            abbrev = kwargs.get("meta_data")[idx].get("Label")
            formula = formula.replace(abbrev, value)
        try:
            result = eval(formula)
        except ZeroDivisionError:
            log.warning(
                f"Formula {formula} for {plate.barcode}{well.hr_position} resulted in a ZeroDivisionError. Setting result to zero."
            )
            result = 0
        except Exception as e:
            log.error(
                f"Error while evaluating formula {formula} for {plate.barcode}{well.hr_position}: {e} "
            )
        return result

    def map(self, data: list[dict], **kwargs) -> None:
        def __debug(msg):
            if kwargs.get("debug"):
                message(msg, "debug", kwargs.get("room_name", None))

        barcode = kwargs.get("barcode")
        try:
            plate = Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:
            if kwargs.get("create_missing_plates"):
                message(
                    f"Plate with barcode {barcode} does not exist. Creating it.",
                    "warning",
                    kwargs.get("room_name", None),
                )

                barcode_specification, _ = BarcodeSpecification.objects.get_or_create(
                    prefix=barcode.split("_")[0],
                    sides=["North"],
                    number_of_plates=4,
                    experiment=Experiment.objects.get(
                        name=kwargs.get("experiment_name")
                    ),
                )
                plate = Plate.objects.create(
                    barcode=barcode,
                    dimension=PlateDimension.by_num_wells(len(data)),
                    experiment=barcode_specification.experiment,
                )
            else:
                message(
                    f"Plate with barcode {barcode} does not exist.",
                    "error",
                    kwargs.get("room_name", None),
                )
                raise ValueError(f"Plate with barcode {barcode} does not exist.")

        with tqdm(
            desc="Processing measurements",
            unit="measurement",
            total=len(data),
        ) as mbar:
            assignment = self.create_measurement_assignment(
                plate,
                kwargs.get("filename"),
            )

            for entry in data:
                __debug(f"Entry: {entry}")
                position = plate.dimension.position(entry.get("position"))
                well = plate.well_at(position)
                if not well:
                    well = Well.objects.create(plate=plate, position=position)

                if not kwargs.get("evaluation"):
                    for idx, value in enumerate(entry.get("values")):
                        label = kwargs.get("meta_data")[idx].get("Label")

                        Measurement.objects.update_or_create(
                            well=well,
                            label=label,
                            measured_at=kwargs.get("measurement_date"),
                            defaults={
                                "value": value,
                                "identifier": entry.get("identifier"),
                                "measurement_assignment": assignment,
                            },
                        )

                else:
                    evaluation_formulas = kwargs.get("evaluation").split(",")
                    measurement_names = kwargs.get("measurement_name").split(",")
                    for idx, evaluation_formula in enumerate(evaluation_formulas):
                        label = measurement_names[idx]
                        result = self.apply_evaluation_formula(
                            evaluation_formula, plate, well, entry, **kwargs
                        )
                        Measurement.objects.update_or_create(
                            well=well,
                            label=label,
                            measured_at=kwargs.get("measurement_date"),
                            defaults={
                                "value": result,
                                "identifier": entry.get("identifier"),
                                "measurement_assignment": assignment,
                            },
                        )
                        plate_mapping = PlateMapping.objects.get(target_plate=plate)
                        plate_mapping.evaluation = evaluation_formula
                        plate_mapping.save()

                mbar.update(1)

    def create_measurement_assignment(self, plate, filename):
        with open(filename, "rb") as file:
            assignment, _ = MeasurementAssignment.objects.update_or_create(
                status="success",
                plate=plate,
                filename=filename,
                measurement_file=File(file, os.path.basename(file.name)),
            )
            return assignment
