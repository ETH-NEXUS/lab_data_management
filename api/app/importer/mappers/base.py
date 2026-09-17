"""
The base class of all mappers.

A mapper reads instrument files (`parse`) and writes their content into the
database (`map`). BaseMapper goes through the files and holds what the mappers
share: creating plates, barcode specifications and measurement assignments.
"""

import os
from contextlib import redirect_stderr
from glob import glob
from typing import Any

from chardet.universaldetector import UniversalDetector
from django.core.management.base import CommandError
from django.core.files import File
from django.core.files.storage import default_storage
from django.db import transaction

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
from importer.command_output import error_text
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

    def __init__(self) -> None:
        # The copies of report files saved in the media folder while one file is
        # mapped, e.g. ["20240610-121212_demo_1_Kx3dP0a.asc"]. A rollback does not
        # delete them, so a file that fails deletes them itself.
        self.stored_files: list[str] = []

    @staticmethod
    def get_files(glob_pattern: str) -> list[str]:
        """All files that match the pattern; "**" also looks into sub folders."""
        return glob(glob_pattern, recursive=True)

    def run(self, glob_pattern: str, **kwargs) -> None:
        """
        Parses and maps every file that matches the pattern, then refreshes
        the materialized views once.

        Every file is mapped in its own transaction: a file with an error
        stores nothing, the error is shown, and the next file is mapped.

        All files share the same kwargs: what one file adds ("xml_file",
        "filename") is still there when the next file is read.
        """
        room_name = kwargs.get("room_name")
        filenames = self.get_files(glob_pattern)
        if not filenames:
            message(f"No files found that match {glob_pattern}.", "warning", room_name)

        failed_files = []
        for filename in filenames:
            message(f"Processing file {filename}...", "info", room_name)
            self.stored_files = []
            try:
                with transaction.atomic():
                    data = self.read_file(filename, kwargs)
                    kwargs.update({"filename": filename})
                    self.map(data, **kwargs)
            except Exception as error:
                self.delete_stored_files()
                failed_files.append(filename)
                message(
                    f"{filename} was not mapped, nothing of it was stored: "
                    f"{error_text(error)}",
                    "error",
                    room_name,
                )
                if not isinstance(error, CommandError):
                    logger.exception(f"Mapping {filename} failed")

        # With a single file, its own error already says everything
        if len(filenames) > 1 and failed_files:
            message(
                f"{len(failed_files)} of {len(filenames)} files were not mapped: "
                f"{', '.join(failed_files)}",
                "error",
                room_name,
            )

        message("Refreshing materialized views...", "info", room_name)
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        ExperimentDetail.refresh(concurrently=True)

    def read_file(self, filename: str, kwargs: dict):
        """
        Parses one file and returns what `parse` returned.

        Changes kwargs in place: sets "xml_file" for files that are opened here.
        """
        if filename.endswith(FILES_PARSED_BY_NAME):
            logger.info(f"Parsing {filename} by its file name")
            return self.parse(filename, **kwargs)

        kwargs.update({"xml_file": filename.endswith(".xml")})
        with open(filename, "r", encoding=detect_encoding(filename)) as file:
            return self.parse(file, **kwargs)

    def parse(self, file: Any, **kwargs) -> Any:
        """Every mapper decides what `file` is and what it returns."""
        raise NotImplementedError

    def map(self, data: Any, **kwargs) -> None:
        """`data` is what `parse` of the same mapper returned."""
        raise NotImplementedError

    def create_measurement_assignment(
        self, plate: Plate, filename: str
    ) -> MeasurementAssignment:
        """Links the measurement file to the plate, with status "success"."""
        with open(filename, "rb") as file:
            assignment, created = MeasurementAssignment.objects.update_or_create(
                status="success",
                plate=plate,
                filename=filename,
                measurement_file=File(file, os.path.basename(file.name)),
            )
        if created:
            # An assignment that was there already keeps the copy it already has
            self.stored_files.append(assignment.measurement_file.name)
        return assignment

    def delete_stored_files(self) -> None:
        """Deletes the report copies of a file whose mapping was rolled back."""
        for name in self.stored_files:
            default_storage.delete(name)
        self.stored_files = []

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
        `experiment_name`; without an experiment name a CommandError is raised.
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
                raise CommandError(text)

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
                plate_name, plate_type, source_plate_name
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
    ) -> PlateDimension:
        """
        The plate dimension for the well count found in the plate names,
        e.g. "Greiner_384PS_781904" -> the dimension with 16 rows and 24 columns.
        """
        names = f"{plate_type} {plate_name} {source_plate_name}"
        try:
            rows, cols = row_col_from_name(names)
            return PlateDimension.objects.get(rows=rows, cols=cols)
        except ValueError as error:
            raise CommandError(
                f"Could not determine plate dimensions for {plate_name}: {error}"
            )
        except PlateDimension.DoesNotExist:
            raise CommandError(f"No plate dimension found for {plate_name}")
