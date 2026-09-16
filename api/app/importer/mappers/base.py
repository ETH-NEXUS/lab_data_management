"""
The base class of all mappers: runs over the files and creates plates.
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
            if filename.endswith(".xlsx") or filename.endswith(".txt"):
                logger.info(f"Processing {filename} as Excel file")
                data = self.parse(filename, **kwargs)
            else:
                xml_file = filename.endswith(".xml")
                # add xml_file to **kwargs
                kwargs.update({"xml_file": xml_file})
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

    def parse(self, file: TextIOWrapper | str | TextIO, **kwargs):
        raise NotImplementedError

    def map(self, data: list[dict], **kwargs) -> None:
        raise NotImplementedError

    def create_measurement_assignment(self, plate, filename):
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
        **kwargs,
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
            dimension=self.get_plate_dimension(
                plate_name, plate_type, source_plate_name, kwargs.get("room_name")
            ),
        )

    def get_plate_dimension(
        self, plate_name: str, plate_type: str, source_plate_name, room_name
    ):
        try:
            rows, cols = row_col_from_name(
                f"{plate_type} {plate_name} {source_plate_name}"
            )
            plate_dimension = PlateDimension.objects.get(rows=rows, cols=cols)
            return plate_dimension
        except ValueError:
            message(
                f"Could not determine plate dimensions for {plate_name}",
                "error",
                room_name,
            )
            raise
        except PlateDimension.DoesNotExist:
            message(f"No plate dimension found", "error", room_name)
            raise ValueError(f"No plate dimension found for {plate_name}")
