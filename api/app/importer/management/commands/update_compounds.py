"""
This script updates the "data" field of compounds in the database using
information from a CSV file provided by the supplier. The CSV file is used
because the data from SDF files is incomplete.
"""

import csv
from typing import Dict
from django.core.management import BaseCommand
from django.core.management.base import CommandError
from django.db import transaction
from compoundlib.models import Compound
from helpers.logger import logger

COMPOUND_NAME_FIELD = "CompoundName"


class Command(BaseCommand):
    help = "Update compound data in the database using a CSV file."

    def add_arguments(self, parser) -> None:
        parser.add_argument(
            "--input_file",
            "-i",
            type=str,
            required=True,
            help="The input file with the compound data, columns separated by tabs.",
        )

    def read_csv_compound_data(self, input_file: str) -> list[dict]:
        """
        Reads compound data from the provided CSV file.

        :param input_file: The path to the input CSV file.
        :return: List of dictionaries containing the CSV data.
        """
        logger.info(f"Reading compound data from {input_file}")
        try:
            # "utf-8-sig" removes the BOM that Excel writes before the first column name
            with open(input_file, "r", encoding="utf-8-sig") as file:
                reader = csv.DictReader(file, delimiter="\t")
                return [row for row in reader]
        except OSError as error:
            raise CommandError(f"Cannot read the file {input_file}: {error}")

    def update_compound_data(
        self, compound: Compound, new_data: Dict[str, str]
    ) -> None:
        """
        Updates the compound's data field with the new data from the CSV.

        :param compound: The compound instance to update.
        :param new_data: Dictionary containing the new data.
        """
        compound_data = compound.data if compound.data else {}
        logger.info(f"Updating compound '{compound.name}' with new data.")
        compound_data.update(new_data)
        compound.data = compound_data
        compound.save()

    def import_compound_data(self, input_file: str) -> None:
        """
        Imports compound data from the CSV and updates the database.

        :param input_file: The path to the input CSV file.
        """
        data = self.read_csv_compound_data(input_file)
        if not data:
            logger.warning("No data found in the CSV file.")
            return

        logger.debug(f"Found CSV headers: {data[0].keys()}")

        for item in data:
            compound_name = item.get(COMPOUND_NAME_FIELD)
            if not compound_name:
                logger.warning(f"Missing '{COMPOUND_NAME_FIELD}' in row: {item}")
                continue

            compound = Compound.objects.filter(name=compound_name).first()
            if compound is None:
                logger.error(f"Compound '{compound_name}' not found in the database.")
                continue

            self.update_compound_data(compound, item)

    def handle(self, *args, **options) -> None:
        """
        Main handler method for the command.
        """
        input_file: str = options.get("input_file")
        logger.info(f"Starting compound data import from {input_file}")
        # All compounds are updated together: after an error nothing is stored
        with transaction.atomic():
            self.import_compound_data(input_file)
        logger.info("Compound data import completed.")
