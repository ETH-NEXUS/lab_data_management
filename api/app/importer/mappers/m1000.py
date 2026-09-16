"""
Maps M1000 reader files (.asc) to measurements.
"""

import os
import re
from datetime import datetime as dt
from io import TextIOWrapper

import pandas as pd
from django.utils import timezone as tz
from tqdm import tqdm

from core.models import (
    BarcodeSpecification,
    Experiment,
    MappingError,
    Measurement,
    Plate,
    PlateDimension,
    Well,
)
from importer.helper import message
from importer.mappers.base import BaseMapper
from importer.mappers.values import convert_sci_to_float


class M1000Mapper(BaseMapper):
    RE_FILENAME = r"(?:(?P<date>[0-9]+)-(?P<time>[0-9]+)_)?(?P<barcode>[^\.]+)\.asc"

    RE_POS = r"^[A-Z]+[0-9]+$"
    RE_ID = r"^[^_]+_[^_]+$"
    RE_NUM = r"^[0-9\.]+$"
    RE_TAB = r"[\t\s]+"
    RE_SCIENTIFIC = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"

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
                f"File name {os.path.basename(file.name)} does not match conventions."
            )

        results = []
        measurement_date = tz.now()
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
                        convert_sci_to_float(part)
                        for idx, part in enumerate(parts)
                        if idx not in [pos_index, id_index]
                        and (
                            re.match(self.RE_NUM, part)
                            or re.match(self.RE_SCIENTIFIC, part)
                        )
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

    def map(self, data: list[dict], **kwargs) -> None:
        def __debug(msg):
            if kwargs.get("debug"):
                message(msg, "debug", kwargs.get("room_name", None))

        barcode = kwargs.get("barcode")
        try:
            plate = Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:

            message(
                f"Plate with barcode {barcode} does not exist. Creating it.",
                "warning",
                kwargs.get("room_name", None),
            )

            barcode_specification, _ = BarcodeSpecification.objects.get_or_create(
                prefix=barcode.split("_")[0],
                sides=["North"],
                number_of_plates=4,
                experiment=Experiment.objects.get(name=kwargs.get("experiment_name")),
            )
            plate = Plate.objects.create(
                barcode=barcode,
                dimension=PlateDimension.by_num_wells(len(data)),
                experiment=barcode_specification.experiment,
            )

        with tqdm(
            desc="Processing measurements",
            unit="measurement",
            total=len(data),
        ) as mbar:
            assignment = self.create_measurement_assignment(
                plate, kwargs.get("filename")
            )

            for entry in data:
                __debug(f"Entry: {entry}")
                position = plate.dimension.position(entry.get("position"))
                well = plate.well_at(position)
                if not well:
                    well = Well.objects.create(plate=plate, position=position)

                measurement_names = None
                if kwargs.get("measurement_name"):
                    measurement_names = kwargs.get("measurement_name").split(",")
                for idx, value in enumerate(entry.get("values")):
                    if measurement_names:
                        label = measurement_names[idx]
                    else:
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
                mbar.update(1)
