"""
Maps Echo transfer reports (CSV or XML) onto plates.
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
        if kwargs.get("xml_file"):
            tree = ET.parse(file)
            root = tree.getroot()
            result = []
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
            for w in root.find("printmap"):
                result.append(
                    {
                        "source_plate_name": source_plate_name,
                        "source_plate_barcode": source_plate_barcode,
                        "destination_plate_name": destination_plate_name,
                        "destination_plate_barcode": destination_plate_barcode,
                        "source_well": w.get("n"),
                        "destination_well": w.get("dn"),
                        "actual_volume": w.get("vl"),
                        "current_fluid_volume": w.get("cvl"),
                        "DMSO": w.get("fc"),
                        "transfer_status": "",
                    }
                )

            return result

        headers = kwargs.get("headers", EchoMapper.DEFAULT_COLUMNS)
        file = self.__fast_forward_to_header_row(file, headers)
        results = []
        # Transfers skipped because some of their required fields are empty.
        # Example: ["Drug08_J A3 -> 2026Wagner12 A3 (Actual Volume empty)"]
        skipped_transfers = []
        reader = csv.DictReader(file, delimiter=",")

        for row in reader:
            # If there are None or empty values in any of the following keys
            # or if there is a second header column we continue.
            # Empty values come from section lines like "[DETAILS],,,,": all
            # columns are present, but only the first one has text.
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
                    row[headers.get(key)] in (None, "")
                    or row[headers.get(key)] == headers.get(key)
                    for key in must_keys
                ]
            ):
                empty_columns = [
                    headers.get(key) for key in must_keys if row[headers.get(key)] == ""
                ]
                # A section line like "[DETAILS],,,," has all required fields
                # empty and is skipped silently. Only a transfer with some empty
                # fields is reported, because it may be a real transfer we lose.
                if 0 < len(empty_columns) < len(must_keys):
                    source = (
                        f"{row[headers.get('source_plate_barcode')]} "
                        f"{row[headers.get('source_well')]}"
                    )
                    destination = (
                        f"{row[headers.get('destination_plate_barcode')]} "
                        f"{row[headers.get('destination_well')]}"
                    )
                    skipped_transfers.append(
                        f"{source} -> {destination} ({', '.join(empty_columns)} empty)"
                    )
                continue

            res_dict = {}
            for new_key, old_key in headers.items():
                if old_key in row:
                    res_dict[new_key] = row[old_key]
            results.append(res_dict)

        # One message for the whole file: the management page only shows the
        # latest message, so one message per row could be overwritten.
        if skipped_transfers:
            message(
                f"Skipped {len(skipped_transfers)} Echo transfers with empty "
                f"required fields: {'; '.join(skipped_transfers)}",
                "warning",
                kwargs.get("room_name", None),
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
                source_plate_name = entry["source_plate_name"]
                source_plate_barcode = entry["source_plate_barcode"]
                destination_plate_name = entry["destination_plate_name"]
                current_fluid_volume = (
                    float(entry["current_fluid_volume"])
                    if entry["current_fluid_volume"]
                    else None
                )
                current_dmso = (
                    float(entry["DMSO"].replace("%", "")) if entry["DMSO"] else None
                )
                if "destination_plate_type" in entry:
                    destination_plate_type = entry["destination_plate_type"]
                else:
                    destination_plate_type = ""
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
                        destination_plate = self.create_plate_by_name_and_barcode(
                            destination_plate_name,
                            destination_plate_type,
                            destination_plate_barcode,
                            source_plate_name,
                            **kwargs,
                        )
                mapping_list_index = (
                    f"{source_plate_barcode}__**__{destination_plate_barcode}"
                )
                if mapping_list_index in mapping_lists:
                    mapping_list = mapping_lists.get(mapping_list_index)
                else:
                    mapping_list = MappingList(target=destination_plate)
                    mapping_lists[mapping_list_index] = mapping_list

                source_well = entry["source_well"]
                destination_well = entry["destination_well"]
                from_pos = source_plate.dimension.position(source_well)
                to_pos = destination_plate.dimension.position(destination_well)

                mapping = Mapping(
                    from_pos=from_pos,
                    to_pos=to_pos,
                    amount=float(entry["actual_volume"]),
                    status=entry["transfer_status"],
                    map_type=source_plate.is_control_plate,  # if it is true, we need to map the type
                    current_amount=current_fluid_volume,
                    current_dmso=current_dmso,
                )
                mapping_list.add(mapping)
                mbar.update(1)

        for mapping_list_index, mapping_list in mapping_lists.items():
            source_plate_barcode = mapping_list_index.split("__**__")[0]
            source_plate = plates.get(source_plate_barcode)
            message(
                f"Mapping {source_plate.barcode} -> {mapping_list.target.barcode}",
                "info",
                kwargs.get("room_name", None),
            )

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
                PlateDetail.refresh(concurrently=True)
                WellDetail.refresh(concurrently=True)
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
