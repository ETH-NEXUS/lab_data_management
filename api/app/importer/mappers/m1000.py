"""
Maps M1000 reader files (.asc) to measurements.

An .asc file has one line per well: the well position, an identifier and one or
more values, e.g. "A1<TAB>SM1_1<TAB>15". The footer at the end has the date of
the measurement, the plate description and the settings of every measured label.
"""

import os
import re
from datetime import datetime as dt
from io import TextIOWrapper
from typing import TypedDict

from django.core.management.base import CommandError
from django.utils import timezone as tz
from tqdm import tqdm

from core.models import Measurement
from importer.helper import message
from importer.mappers.base import BaseMapper
from importer.mappers.values import convert_sci_to_float

# In the footer, the key that names a measured label
META_DATA_LABEL = "Label"


class M1000Entry(TypedDict):
    """One value line of an .asc file."""

    position: str  # the well, e.g. "A1"
    identifier: str  # e.g. "SM1_1"
    values: list[float | None]  # one value per measured label, e.g. [15.0]


class M1000Data(TypedDict):
    """What `parse` returns for one .asc file."""

    barcode: str  # from the file name, e.g. "demo_1"
    measurement_date: dt  # e.g. datetime(2011, 11, 11, 11, 11, 11)
    plate_description: str | None  # e.g. "test"
    # The settings of every label, in the order of the value columns,
    # e.g. [{"Label": "Label1", "Integration time": "1000 ms"}]
    meta_data: list[dict[str, str]]
    entries: list[M1000Entry]


def debug_message(text: str, kwargs: dict) -> None:
    """Sends a debug message, but only when the mapper runs with debug=True."""
    if kwargs.get("debug"):
        message(text, "debug", kwargs.get("room_name"))


class M1000Mapper(BaseMapper):
    # File name with optional date and time, e.g. "20240610-121212_demo_1.asc"
    RE_FILENAME = r"(?:(?P<date>[0-9]+)-(?P<time>[0-9]+)_)?(?P<barcode>[^\.]+)\.asc"

    # A well position, e.g. "A1" or "AB12"
    RE_POS = r"^[A-Z]+[0-9]+$"
    # An identifier with exactly one underscore, e.g. "SM1_1"
    RE_ID = r"^[^_]+_[^_]+$"
    # A plain number, e.g. "15" or "1.5"
    RE_NUM = r"^[0-9\.]+$"
    # What separates the columns: tabs and spaces
    RE_TAB = r"[\t\s]+"
    # A number, also in scientific notation, e.g. "1.5E+02". Used with re.match,
    # so only the start of the text has to look like a number.
    RE_SCIENTIFIC = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"

    # Footer: "Date of measurement: 2011-11-11/Time of measurement: 11:11:11"
    RE_DATE_OF_MEASUREMENT = r"^Date of measurement: (?P<date>[^\/]+)\/Time of measurement: (?P<time>[0-9:]+)$"
    # Footer: "Plate Description: test"
    RE_PLATE_DESCRIPTION = r"^Plate Description: (?P<description>.+)$"
    # Footer: a label setting, indented by four spaces, e.g. "    Label: Label1"
    RE_META_DATA = r"^    (?P<key>[^:]+): (?P<value>.+)$"

    def determine_indexes(self, file: TextIOWrapper) -> tuple[int, int]:
        """
        Finds the column of the position and the column of the identifier; the
        other columns are values. Returns e.g. (0, 1) for "A1<TAB>SM1_1<TAB>15".

        The first line with exactly one column that looks like a position and
        exactly one column that looks like an identifier decides. Afterwards
        the file is read again from its start.
        """
        for line in file:
            parts = re.split(self.RE_TAB, line.strip())
            position_columns = []
            identifier_columns = []
            for index, part in enumerate(parts):
                if re.match(self.RE_POS, part):
                    position_columns.append(index)
                if re.match(self.RE_ID, part):
                    identifier_columns.append(index)

            if len(position_columns) == 1 and len(identifier_columns) == 1:
                file.seek(0)
                return position_columns[0], identifier_columns[0]

        file.seek(0)
        raise CommandError(f"File has not the desired format: {file.name}")

    def parse(self, file: TextIOWrapper, **kwargs) -> M1000Data:
        """
        Reads the value lines and the footer.

        Returned data example:
        {"barcode": "demo_1", "plate_description": "test",
         "measurement_date": datetime(2011, 11, 11, 11, 11, 11),
         "meta_data": [{"Label": "Label1", "Integration time": "1000 ms"}],
         "entries": [{"position": "A1", "identifier": "SM1_1", "values": [15.0]}]}
        """
        barcode = self.barcode_from_file_name(file.name)

        entries = []
        # Used when the footer has no date of measurement
        measurement_date = tz.now()
        plate_description = None
        # The settings of every label, one dict per label. A key that appears
        # again starts the settings of the next label.
        meta_data: list[dict[str, str]] = [{}]
        position_column, identifier_column = self.determine_indexes(file)

        for line in file:
            parts = re.split(self.RE_TAB, line)
            # A value line has at least three parts (position, identifier, value)
            if len(parts) >= 3 and re.match(self.RE_POS, parts[position_column]):
                entry = self.read_value_line(parts, position_column, identifier_column)
                # A line without values is ignored
                if len(entry["values"]) > 0:
                    debug_message(f"result: {entry}", kwargs)
                    entries.append(entry)
            elif match := re.match(self.RE_DATE_OF_MEASUREMENT, line):
                measurement_date = dt.strptime(
                    f"{match.group('date')} {match.group('time')}", "%Y-%m-%d %H:%M:%S"
                )
            elif match := re.match(self.RE_PLATE_DESCRIPTION, line):
                plate_description = match.group("description")
            elif match := re.match(self.RE_META_DATA, line):
                key = match.group("key")
                if key in meta_data[-1]:
                    meta_data.append({})
                meta_data[-1][key] = match.group("value")

        return {
            "barcode": barcode,
            "measurement_date": measurement_date,
            "plate_description": plate_description,
            # The labels are listed in the footer in the reverse order of the value columns
            "meta_data": list(reversed(meta_data)),
            "entries": entries,
        }

    def barcode_from_file_name(self, path: str) -> str:
        """ "/data/20240610-121212_demo_1.asc" -> "demo_1"."""
        file_name = os.path.basename(path)
        match = re.match(self.RE_FILENAME, file_name)
        if not match:
            raise CommandError(f"File name {file_name} does not match conventions.")
        return match.group("barcode")

    def read_value_line(
        self, parts: list[str], position_column: int, identifier_column: int
    ) -> M1000Entry:
        """
        ["A1", "SM1_1", "15", ""] -> {"position": "A1", "identifier": "SM1_1", "values": [15.0]}.
        Every other column that looks like a number becomes a value. A column
        that only starts like a number (e.g. "12abc") refuses the file, because
        a measurement needs a value.
        """
        values = []
        for index, part in enumerate(parts):
            if index in (position_column, identifier_column):
                continue
            if re.match(self.RE_NUM, part) or re.match(self.RE_SCIENTIFIC, part):
                value = convert_sci_to_float(part)
                if value is None:
                    raise CommandError(
                        f"The value '{part}' of well {parts[position_column]} is not "
                        "a number."
                    )
                values.append(value)

        return {
            "position": parts[position_column],
            "identifier": parts[identifier_column],
            "values": values,
        }

    def map(self, data: M1000Data, **kwargs) -> None:
        """
        Stores every value of every entry as a measurement of its well, and
        links the file to the plate with a measurement assignment.
        """
        entries = data["entries"]
        labels = self.measurement_labels(
            entries, data["meta_data"], kwargs.get("measurement_name")
        )
        plate = self.find_or_create_measured_plate(
            data["barcode"],
            len(entries),
            kwargs.get("room_name"),
            kwargs.get("experiment_name"),
        )

        with tqdm(
            desc="Processing measurements",
            unit="measurement",
            total=len(entries),
        ) as progress:
            assignment = self.create_measurement_assignment(plate, kwargs["filename"])

            for entry in entries:
                debug_message(f"Entry: {entry}", kwargs)
                position = plate.dimension.position(entry.get("position"))
                well = plate.well_at(position, create_if_not_exist=True)

                for index, value in enumerate(entry.get("values")):
                    Measurement.objects.update_or_create(
                        well=well,
                        label=labels[index],
                        measured_at=data["measurement_date"],
                        defaults={
                            "value": value,
                            "identifier": entry.get("identifier"),
                            "measurement_assignment": assignment,
                        },
                    )
                progress.update(1)

    @staticmethod
    def measurement_labels(
        entries: list[M1000Entry],
        meta_data: list[dict[str, str]],
        measurement_name: str | None,
    ) -> list[str]:
        """
        The label of every value column: the names of `measurement_name`
        ("Lum,Fluo" -> ["Lum", "Fluo"]), or else the "Label" of the settings in
        the footer of the file. A file without a label for every column is refused.
        """
        columns = max((len(entry["values"]) for entry in entries), default=0)
        if measurement_name:
            labels = measurement_name.split(",")
            if len(labels) < columns:
                raise CommandError(
                    f"The file has {columns} value columns, but only "
                    f"{len(labels)} measurement names were given: {measurement_name}."
                )
            return labels

        labels = [settings.get(META_DATA_LABEL, "") for settings in meta_data]
        if len(labels) < columns or not all(labels[:columns]):
            raise CommandError(
                f"The file has {columns} value columns, but its footer does not "
                "name every label. Please give the measurement names, "
                "for example 'Lum,Fluo'."
            )
        return labels
