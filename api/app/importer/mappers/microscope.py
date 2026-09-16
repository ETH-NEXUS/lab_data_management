"""
Maps C10 reader and imager files (.xlsx or .txt) to measurements.
"""

import os
import re

from openpyxl import load_workbook
from tqdm import tqdm

from core.models import (
    BarcodeSpecification,
    Experiment,
    Measurement,
    Plate,
    PlateDimension,
    Well,
    WellType,
)
from helpers.logger import logger
from importer.helper import message
from importer.mappers.base import BaseMapper
from importer.mappers.values import convert_sci_to_float, convert_string_to_datetime


class MicroscopeMapper(BaseMapper):
    RE_FILENAME = (
        r"(?P<date>\d+)[-_](?P<time>\d+)[-_](?P<barcode>[^\.]+)\.(?P<ext>xlsx|txt)$"
    )

    def parse(self, file, **kwargs):
        filename = file
        basename = os.path.basename(filename)
        match = re.match(self.RE_FILENAME, basename)
        if match:
            date = match.group("date")
            time = match.group("time")
            barcode = match.group("barcode")
            ext = match.group("ext")
        else:
            message(
                f"Filename {basename} does not match expected pattern.",
                "error",
                kwargs.get("room_name", None),
            )
            raise ValueError(f"Filename {basename} does not match expected pattern.")

        kwargs.update(
            {"filename": filename, "barcode": barcode, "date": date, "time": time}
        )
        measurement_name = kwargs.get("measurement_name", "Lum")

        if ext == "xlsx":
            wb = load_workbook(file)
            sheet = wb.active
            metadata = self.__parse_metadata(sheet)
            results = self.__parse_results(sheet)
            layout = self.__parse_layout(sheet, len(results))
        elif ext == "txt":
            content = open(file, "r")
            lines = [line.strip() for line in content.readlines() if line.strip()]
            metadata = self.__parse_metadata_txt(lines)
            results = self.__parse_results_txt(lines, measurement_name)
            layout = (
                {}
            )  # Implement __parse_layout_txt if layout info is present in .txt files
            # Prefer date and time from metadata if available
            date = metadata.get("Date", date)
            time = metadata.get("Time", time)
            logger.info(f"Date: {date}, Time: {time}")
            content.close()
        else:
            raise ValueError(f"Unsupported file extension: {ext}")

        return {
            "metadata": metadata,
            "results": results,
            "date": date,
            "time": time,
            "barcode": barcode,
            "layout": layout,
        }

    def map(self, data: dict, **kwargs) -> None:

        RE_NUMBER = r"^[0-9]+(\.[0-9]+)?$"
        RE_SCIENCE = r"^[0-9\.]+[eE][+-]?[0-9]+$"

        barcode = data["barcode"]
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
                dimension=PlateDimension.by_num_wells(len(data["results"])),
                experiment=barcode_specification.experiment,
            )
        measured_at = convert_string_to_datetime(data["date"], data["time"])
        with tqdm(
            desc="Processing microscope output",
            unit="measurement",
            total=len(data["results"]),
        ) as mbar:
            for entry in data["results"]:
                if not entry.get("Well") or entry.get("Well") == "Well":
                    continue

                position = plate.dimension.position(entry.get("Well"))

                well = plate.well_at(position)
                if not well:

                    well = Well.objects.create(plate=plate, position=position)
                if data["layout"] and entry.get("Well") in data["layout"]:
                    well_type = data["layout"][entry.get("Well")]
                    well.type = WellType.objects.get(name=well_type)
                    well.save()
                for key, value in entry.items():
                    if key in ["Well ID", "Well"]:
                        continue
                    if re.match(RE_NUMBER, str(value)):
                        value = float(value)
                    elif re.match(RE_SCIENCE, str(value)):
                        value = convert_sci_to_float(value)
                    else:
                        continue
                    Measurement.objects.update_or_create(
                        well=well,
                        label=key,
                        measured_at=measured_at,
                        defaults={
                            "value": value,
                        },
                    )

                mbar.update(1)
        self.create_measurement_assignment(plate, kwargs.get("filename"))

    def __parse_metadata(self, sheet):
        metadata = {}
        current_label = None
        for row in sheet.iter_rows(min_row=1, max_row=38):
            for index, cell in enumerate(row):
                if index == 0 and cell.value:
                    current_label = str(cell.value).strip().lstrip().replace(":", "")
                    metadata[current_label] = []

                elif index != 0 and cell.value:
                    str_value = str(cell.value)
                    if ": " in str_value:
                        k = str_value.split(":")[0].strip().lstrip()
                        v = str_value.split(":")[1].strip().lstrip()
                        metadata[k] = v
                    else:
                        metadata[current_label].append(str_value)
                if str(cell.value).strip().lstrip().lower() == "results":
                    break

        return metadata

    def __parse_results(self, sheet):
        results_data = []
        headers = []
        results_start_row = None
        for index, row in enumerate(sheet.iter_rows(values_only=True)):
            if (
                row[1]
                and str(row[1]).lower() == "well id"
                or str(row[1]).lower() == "well"
            ):
                for cell in sheet.iter_rows(
                    min_row=index + 1,
                    max_row=index + 1,
                    values_only=True,  # we use index +1 because indices start with 1 in Excel
                ):
                    if cell:
                        headers.append(cell)
                results_start_row = index + 2
                break
        headers = [header for header in headers[0] if header is not None]
        if headers and results_start_row:
            for row in sheet.iter_rows(min_row=results_start_row, values_only=True):
                if not any(row):
                    continue
                row_data = dict(zip(headers, row[1:]))
                results_data.append(row_data)

        return results_data

    def __parse_layout(self, sheet, plate):
        layout_start_row = None
        layout_end_row = None
        layout_data = []
        position_type = {}
        for index, row in enumerate(sheet.iter_rows(values_only=True)):
            if row[0] and str(row[0]).lower() == "layout":
                layout_start_row = index + 2
            if row[0] and str(row[0]).lower() == "results":
                layout_end_row = index - 1
        if layout_start_row and layout_end_row:
            for row in sheet.iter_rows(
                min_row=layout_start_row, max_row=layout_end_row, values_only=True
            ):
                if not any(row):
                    continue
                layout_data.append(row[1:])

        for item in layout_data:
            if item[0]:
                for index, value in enumerate(item):
                    if value in ["POS", "NEG"]:
                        well_position = f"{item[0]}{index}"
                        well_type = "P" if value == "POS" else "N"
                        position_type[well_position] = well_type

        return position_type

    def __parse_metadata_txt(self, lines):
        metadata = {}
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line == "Results":
                break  # Stop parsing metadata when 'Results' is reached
            if not line:
                i += 1
                continue
            if "\t" in line:
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    key, value = parts
                    metadata[key.strip()] = value.strip()
            elif ":" in line:
                parts = line.split(":", 1)
                if len(parts) == 2:
                    key, value = parts
                    metadata[key.strip()] = value.strip()

            i += 1

        return metadata

    def __parse_results_txt(self, lines, measurement_name="Lum"):
        results = []
        # Find the index where 'Results' section starts
        i = 0
        while i < len(lines):
            if lines[i] == "Results":
                i += 1  # Skip the 'Results' header
                break
            i += 1

        # Now parse the results
        while i < len(lines):
            line = lines[i]
            if line.strip() == "":
                break
            parts = line.split("\t")
            if len(parts) == 2:
                well, lum = parts
                if well.strip() == "Well":
                    i += 1
                    continue
                results.append({"Well": well.strip(), measurement_name: lum.strip()})
            i += 1

        return results
