"""
The base class of all mappers.

A mapper reads instrument files (`parse`) and writes their content into the
database (`map`). BaseMapper goes through the files and holds what the mappers
share: creating plates, barcode specifications and measurement assignments.
"""

import os
from contextlib import redirect_stderr
from glob import glob
from io import TextIOWrapper
from typing import TextIO

from chardet.universaldetector import UniversalDetector
from django.core.files import File

from core.models import (
    BarcodeSpecification,
    Experiment,
    ExperimentDetail,
    MeasurementAssignment,
    Plate,
    PlateDetail,
    PlateDimension,
    WellDetail,
)
from helpers.logger import logger
from importer.helper import message, row_col_from_name

# A barcode specification that the importer creates on its own gets these values.
NEW_BARCODE_SPECIFICATION_SIDES = ("North",)
NEW_BARCODE_SPECIFICATION_NUMBER_OF_PLATES = 4

# The mappers read these files themselves (with openpyxl, or line by line), so
# `parse` gets the file name instead of an open file.
FILES_PARSED_BY_NAME = (".xlsx", ".txt")


def detect_encoding(filename: str) -> str | None:
    """Guesses the text encoding of a file, e.g. "ascii" or "utf-8"."""
    # chardet can print warnings about the file content, they are not needed
    with redirect_stderr(None):
        detector = UniversalDetector()
        with open(filename, "rb") as file:
            for line in file:
                detector.feed(line)
                if detector.done:
                    break
            detector.close()
        return detector.result.get("encoding")


def barcode_prefix(barcode: str) -> str:
    """The part of a barcode before the first underscore: "ABC_1" -> "ABC"."""
    return barcode.split("_")[0]


class BaseMapper:
    """
    Subclasses implement:
    - parse(file, **kwargs): reads one file, e.g. into a list of dicts
    - map(data, **kwargs): writes what parse returned into the database
    """

    @staticmethod
    def get_files(glob_pattern: str) -> list[str]:
        """All files that match the pattern; "**" also looks into sub folders."""
        return glob(glob_pattern, recursive=True)

    def run(self, glob_pattern: str, **kwargs) -> None:
        """
        Parses and maps every file that matches the pattern, then refreshes
        the materialized views once.

        All files share the same kwargs: what one file adds ("xml_file",
        "filename", the extra information of an M1000 file) is still there
        when the next file is read.
        """
        for filename in self.get_files(glob_pattern):
            message(f"Processing file {filename}...", "info", kwargs.get("room_name"))
            data = self.read_file(filename, kwargs)
            kwargs.update({"filename": filename})
            self.map(data, **kwargs)

        message("Refreshing materialized views...", "info", kwargs.get("room_name"))
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

    def read_file(self, filename: str, kwargs: dict):
        """
        Parses one file and returns what `parse` returned.

        Changes kwargs in place: sets "xml_file" for files that are opened
        here, and adds the extra information that `parse` may return.
        """
        encoding = detect_encoding(filename)

        if filename.endswith(FILES_PARSED_BY_NAME):
            logger.info(f"Processing {filename} as Excel file")
            return self.parse(filename, **kwargs)

        kwargs.update({"xml_file": filename.endswith(".xml")})
        with open(filename, "r", encoding=encoding) as file:
            parsed = self.parse(file, **kwargs)

        # The M1000 mapper returns a tuple: (data, extra information)
        if isinstance(parsed, tuple):
            data = parsed[0]
            extra_information = parsed[1]
            kwargs.update(extra_information)
            return data
        return parsed

    def parse(self, file: TextIOWrapper | str | TextIO, **kwargs):
        raise NotImplementedError

    def map(self, data: list[dict], **kwargs) -> None:
        raise NotImplementedError

    def create_measurement_assignment(
        self, plate: Plate, filename: str
    ) -> MeasurementAssignment:
        """Links the measurement file to the plate, with status "success"."""
        with open(filename, "rb") as file:
            assignment, _ = MeasurementAssignment.objects.update_or_create(
                status="success",
                plate=plate,
                filename=filename,
                measurement_file=File(file, os.path.basename(file.name)),
            )
            return assignment

    def create_plate_by_name_and_barcode(
        self,
        plate_name: str,
        plate_type: str,
        barcode: str,
        source_plate_name: str,
        room_name: str | None = None,
        experiment_name: str | None = None,
    ) -> Plate:
        """
        Creates a plate that is not in the database yet.

        The experiment of the plate comes from the barcode specification of the
        barcode prefix. A missing specification is created for the experiment
        `experiment_name`; without an experiment name a ValueError is raised.
        """
        try:
            barcode_specification = BarcodeSpecification.objects.get(
                prefix=barcode_prefix(barcode)
            )
        except BarcodeSpecification.DoesNotExist:
            if not experiment_name:
                text = (
                    f"No barcode specification found for {barcode} and no experiment "
                    "name is provided. Please provide the experiment name in order to "
                    "create the missing barcode specifications."
                )
                message(text, "error", room_name)
                raise ValueError(text)

            message(
                f"No barcode specification found for {barcode}. Creating it.",
                "warning",
                room_name,
            )
            barcode_specification = self.get_or_create_barcode_specification(
                barcode, experiment_name
            )

        return Plate.objects.create(
            barcode=barcode,
            experiment=barcode_specification.experiment,
            dimension=self.get_plate_dimension(
                plate_name, plate_type, source_plate_name, room_name
            ),
        )

    def get_or_create_barcode_specification(
        self, barcode: str, experiment_name: str | None
    ) -> BarcodeSpecification:
        """
        The barcode specification of the barcode prefix in the experiment, with
        the values the importer uses for new specifications.
        """
        barcode_specification, _ = BarcodeSpecification.objects.get_or_create(
            prefix=barcode_prefix(barcode),
            # A new list every time, so the shared default can not be changed
            sides=list(NEW_BARCODE_SPECIFICATION_SIDES),
            number_of_plates=NEW_BARCODE_SPECIFICATION_NUMBER_OF_PLATES,
            experiment=Experiment.objects.get(name=experiment_name),
        )
        return barcode_specification

    def find_or_create_measured_plate(
        self,
        barcode: str,
        number_of_wells: int,
        room_name: str | None,
        experiment_name: str | None,
    ) -> Plate:
        """
        The plate of a measurement file. A missing plate is created for the
        experiment `experiment_name`, with the dimension that fits the number
        of wells in the file (e.g. 384 values -> the 384 well dimension).
        """
        try:
            return Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:
            message(
                f"Plate with barcode {barcode} does not exist. Creating it.",
                "warning",
                room_name,
            )
            barcode_specification = self.get_or_create_barcode_specification(
                barcode, experiment_name
            )
            return Plate.objects.create(
                barcode=barcode,
                dimension=PlateDimension.by_num_wells(number_of_wells),
                experiment=barcode_specification.experiment,
            )

    def get_plate_dimension(
        self,
        plate_name: str,
        plate_type: str,
        source_plate_name: str,
        room_name: str | None,
    ) -> PlateDimension:
        """
        The plate dimension for the well count found in the plate names,
        e.g. "Greiner_384PS_781904" -> the dimension with 16 rows and 24 columns.
        """
        names = f"{plate_type} {plate_name} {source_plate_name}"
        try:
            rows, cols = row_col_from_name(names)
            return PlateDimension.objects.get(rows=rows, cols=cols)
        except ValueError:
            message(
                f"Could not determine plate dimensions for {plate_name}",
                "error",
                room_name,
            )
            raise
        except PlateDimension.DoesNotExist:
            message("No plate dimension found", "error", room_name)
            raise ValueError(f"No plate dimension found for {plate_name}")
