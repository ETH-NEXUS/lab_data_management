"""
Maps Echo transfer reports onto plates.

An Echo report lists liquid transfers, each from a well of a source plate to a
well of a destination plate. `parse` reads the report (CSV or XML) into a list
of transfers. `map` groups the transfers by plate pair and lets Plate.map move
the compounds.
"""

import csv
import os
import xml.etree.ElementTree as ET
from io import TextIOWrapper
from itertools import dropwhile

from django.core.files import File
from tqdm import tqdm

from core.models import Plate, PlateDetail, PlateMapping, WellDetail
from core.utils.plates.mapping import Mapping, MappingList
from importer.config import Config
from importer.helper import message
from importer.mappers.base import BaseMapper

# A row of a CSV report is only a transfer if all these columns have a value.
REQUIRED_COLUMNS = (
    "source_plate_barcode",
    "source_plate_type",
    "source_well",
    "destination_plate_name",
    "destination_plate_barcode",
    "destination_well",
    "actual_volume",
)

# Transfers whose source plate is not in the database are tried again this
# many times, at the end of `map`.
MAX_QUEUE_RETRIES = 3


def number_or_none(text: str) -> float | None:
    """ "10.184" -> 10.184, and "" or None -> None."""
    return float(text) if text else None


def percent_or_none(text: str) -> float | None:
    """ "98.7%" or "98.7" -> 98.7, and "" or None -> None."""
    return float(text.replace("%", "")) if text else None


class EchoMapper(BaseMapper):
    DEFAULT_COLUMNS = Config.current.importer.echo.default.columns

    def parse(self, file: TextIOWrapper, **kwargs) -> list[dict]:
        """
        Reads an Echo report into a list of transfers.

        Returned data example:
        [{"source_plate_name": "384LDV_DMSO", "source_plate_barcode": "Drug08_J",
          "source_well": "L11", "destination_plate_name": "Greiner_384PS_781904",
          "destination_plate_barcode": "2026Wagner12", "destination_well": "L11",
          "actual_volume": "10", "transfer_status": "",
          "current_fluid_volume": "10.184", "DMSO": "98.744"}]
        """
        if kwargs.get("xml_file"):
            return self.parse_xml(file)

        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        return self.parse_csv(file, headers, kwargs.get("room_name"))

    def parse_xml(self, file: TextIOWrapper) -> list[dict]:
        """
        An XML report names the two plates once in <plateInfo>, and every
        transfer is a <w> element in <printmap>.
        """
        root = ET.parse(file).getroot()

        source_plate_name = None
        source_plate_barcode = None
        destination_plate_name = None
        destination_plate_barcode = None
        for plate in root.find("plateInfo"):
            if plate.get("type") == "source":
                source_plate_name = plate.get("name")
                source_plate_barcode = plate.get("barcode")
            elif plate.get("type") == "destination":
                destination_plate_name = plate.get("name")
                destination_plate_barcode = plate.get("barcode")

        transfers = []
        for well in root.find("printmap"):
            transfers.append(
                {
                    "source_plate_name": source_plate_name,
                    "source_plate_barcode": source_plate_barcode,
                    "destination_plate_name": destination_plate_name,
                    "destination_plate_barcode": destination_plate_barcode,
                    "source_well": well.get("n"),
                    "destination_well": well.get("dn"),
                    "actual_volume": well.get("vl"),
                    "current_fluid_volume": well.get("cvl"),
                    "DMSO": well.get("fc"),
                    "transfer_status": "",
                }
            )
        return transfers

    def parse_csv(self, file: TextIOWrapper, headers: dict, room_name) -> list[dict]:
        """
        headers maps our keys to the column names of the report,
        e.g. {"source_well": "Source Well", "DMSO": "% DMSO", ...}.
        """
        rows = csv.DictReader(self.skip_to_header_row(file, headers), delimiter=",")
        transfers = []
        # Example: ["Drug08_J A3 -> 2026Wagner12 A3 (Actual Volume empty)"]
        skipped_transfers = []

        for row in rows:
            if self.is_transfer_row(row, headers):
                transfers.append(self.rename_columns(row, headers))
                continue

            empty_columns = [
                headers.get(key)
                for key in REQUIRED_COLUMNS
                if row[headers.get(key)] == ""
            ]
            # A section line like "[DETAILS],,,," has all required columns
            # empty and is skipped silently. Only a row with some empty columns
            # is reported, because it may be a real transfer that gets lost.
            if 0 < len(empty_columns) < len(REQUIRED_COLUMNS):
                skipped_transfers.append(
                    self.describe_skipped_transfer(row, headers, empty_columns)
                )

        # One message for the whole file: the management page only shows the
        # latest message, so one message per row could be overwritten.
        if skipped_transfers:
            message(
                f"Skipped {len(skipped_transfers)} Echo transfers with empty "
                f"required fields: {'; '.join(skipped_transfers)}",
                "warning",
                room_name,
            )
        return transfers

    @staticmethod
    def skip_to_header_row(file: TextIOWrapper, headers: dict):
        """
        The lines of the file from the header row on. The header row is the
        first line that contains the first column name, e.g. "Source Plate Barcode".
        """
        first_column_name = list(headers.values())[0]
        return dropwhile(lambda line: first_column_name not in line, file)

    @staticmethod
    def is_transfer_row(row: dict, headers: dict) -> bool:
        """
        False for rows that are not transfers: an empty required column (for
        example a section line like "[DETAILS],,,,") or a repeated header row.
        """
        for key in REQUIRED_COLUMNS:
            column_name = headers.get(key)
            value = row[column_name]
            if value in (None, "") or value == column_name:
                return False
        return True

    @staticmethod
    def rename_columns(row: dict, headers: dict) -> dict:
        """
        {"Source Well": "A3", ...} -> {"source_well": "A3", ...}.
        Columns that the report does not have are left out.
        """
        transfer = {}
        for key, column_name in headers.items():
            if column_name in row:
                transfer[key] = row[column_name]
        return transfer

    @staticmethod
    def describe_skipped_transfer(row: dict, headers: dict, empty_columns: list) -> str:
        """Example: "Drug08_J A3 -> 2026Wagner12 A3 (Actual Volume empty)"."""
        source = (
            f"{row[headers.get('source_plate_barcode')]} "
            f"{row[headers.get('source_well')]}"
        )
        destination = (
            f"{row[headers.get('destination_plate_barcode')]} "
            f"{row[headers.get('destination_well')]}"
        )
        return f"{source} -> {destination} ({', '.join(empty_columns)} empty)"

    def map(self, data: list[dict], **kwargs) -> None:
        """
        Groups the transfers by (source plate, destination plate), maps every
        group with Plate.map and stores the report file with the plate mapping.

        Transfers whose source plate does not exist are put into a queue and
        tried again at the end, at most MAX_QUEUE_RETRIES times.
        """
        room_name = kwargs.get("room_name")
        # Source plates by barcode, so every plate is loaded only once
        plates = {}
        # (source plate, MappingList) by (source barcode, destination barcode),
        # in the order the plate pairs appear in the report
        plate_pairs = {}
        # Transfers whose source plate was not found
        queue = []

        with tqdm(
            desc="Processing mappings",
            unit="mappings",
            total=len(data),
        ) as progress:
            for transfer in data:
                # All values are read first: a broken number stops the import
                # here, even for a transfer whose source plate is missing.
                source_plate_name = transfer["source_plate_name"]
                source_plate_barcode = transfer["source_plate_barcode"]
                destination_plate_name = transfer["destination_plate_name"]
                current_amount = number_or_none(transfer["current_fluid_volume"])
                current_dmso = percent_or_none(transfer["DMSO"])
                destination_plate_type = transfer.get("destination_plate_type", "")
                destination_plate_barcode = transfer["destination_plate_barcode"]

                source_plate = self.find_source_plate(
                    source_plate_barcode, plates, room_name
                )
                if source_plate is None:
                    queue.append(transfer)
                    continue

                destination_plate = self.find_or_create_destination_plate(
                    destination_plate_barcode,
                    destination_plate_name,
                    destination_plate_type,
                    source_plate_name,
                    plates,
                    room_name,
                    kwargs.get("experiment_name"),
                )

                pair = (source_plate_barcode, destination_plate_barcode)
                if pair not in plate_pairs:
                    plate_pairs[pair] = (
                        source_plate,
                        MappingList(target=destination_plate),
                    )
                _, mapping_list = plate_pairs[pair]
                mapping_list.add(
                    self.build_mapping(
                        transfer,
                        source_plate,
                        destination_plate,
                        current_amount,
                        current_dmso,
                    )
                )
                progress.update(1)

        self.map_plate_pairs(plate_pairs, kwargs)

        retries = kwargs.get("try_queue", 0)
        if retries < MAX_QUEUE_RETRIES and len(queue) > 0:
            kwargs.update({"try_queue": retries + 1})
            self.map(queue, **kwargs)

    def find_source_plate(self, barcode: str, plates: dict, room_name) -> Plate | None:
        """
        The source plate from the cache or the database, or None (with a
        warning) if it does not exist.
        """
        if barcode in plates:
            return plates.get(barcode)
        try:
            plate = Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:
            # The line break and the spaces belong to the message text: this
            # is how users have always seen it on the management page.
            message(
                f"""Source plate with barcode {barcode} does not exist.
                            I try again later...""",
                "warning",
                room_name,
            )
            return None
        plates[barcode] = plate
        return plate

    def find_or_create_destination_plate(
        self,
        barcode: str,
        plate_name: str,
        plate_type: str,
        source_plate_name: str,
        plates: dict,
        room_name,
        experiment_name: str,
    ) -> Plate:
        """
        The destination plate from the cache or the database. A missing plate
        is created. Destination plates are not added to the cache.
        """
        if barcode in plates:
            return plates.get(barcode)
        try:
            return Plate.objects.get(barcode=barcode)
        except Plate.DoesNotExist:
            message(f"Creating destination plate {plate_name}, {plate_type}")
            return self.create_plate_by_name_and_barcode(
                plate_name,
                plate_type,
                barcode,
                source_plate_name,
                room_name=room_name,
                experiment_name=experiment_name,
            )

    @staticmethod
    def build_mapping(
        transfer: dict,
        source_plate: Plate,
        destination_plate: Plate,
        current_amount: float | None,
        current_dmso: float | None,
    ) -> Mapping:
        """One well-to-well transfer, with the fill level the Echo reported."""
        return Mapping(
            from_pos=source_plate.dimension.position(transfer["source_well"]),
            to_pos=destination_plate.dimension.position(transfer["destination_well"]),
            amount=float(transfer["actual_volume"]),
            status=transfer["transfer_status"],
            # A control plate passes its well types on to the destination plate
            map_type=source_plate.is_control_plate,
            current_amount=current_amount,
            current_dmso=current_dmso,
        )

    def map_plate_pairs(self, plate_pairs: dict, kwargs: dict) -> None:
        """
        Maps every plate pair. A successful mapping is stored as PlateMapping
        together with the report file, and the views are refreshed.
        """
        room_name = kwargs.get("room_name")
        for source_plate, mapping_list in plate_pairs.values():
            target_plate = mapping_list.target
            pair_text = f"{source_plate.barcode} -> {target_plate.barcode}"
            message(f"Mapping {pair_text}", "info", room_name)

            if not source_plate.map(mapping_list, target_plate):
                message(f"Error mapping {pair_text}", "error", room_name)
                continue

            with open(kwargs["filename"], "rb") as file:
                PlateMapping.objects.create(
                    source_plate=source_plate,
                    target_plate=target_plate,
                    mapping_file=File(file, os.path.basename(file.name)),
                )
            message(f"Mapped {pair_text}", "info", room_name)
            PlateDetail.refresh(concurrently=True)
            WellDetail.refresh(concurrently=True)
