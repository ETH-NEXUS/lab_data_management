"""
This script updates the "data" field of compounds in the database using
information from a CSV file provided by the supplier. The CSV file is used
because the data from SDF files is incomplete.
"""
import csv
from typing import Dict, Any
from django.core.management import BaseCommand
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
            help="The input CSV file containing compound data.",
        )

    def read_csv_compound_data(self, input_file: str) -> list[str] | list[Any]:
        """
        Reads compound data from the provided CSV file.

        :param input_file: The path to the input CSV file.
        :return: List of dictionaries containing the CSV data.
        """
        logger.info(f"Reading compound data from {input_file}")
        try:
            with open(input_file, "r") as file:
                reader = csv.DictReader(file, delimiter="\t")
                return [row for row in reader]
        except Exception as e:
            logger.error(f"Failed to read the CSV file {input_file}: {str(e)}")
            return []

    def update_compound_data(
        self, compound: Compound, new_data: Dict[str, str]
    ) -> None:
        """
        Updates the compound's data field with the new data from the CSV.

        :param compound: The compound instance to update.
        :param new_data: Dictionary containing the new data.
        """
        try:
            compound_data = compound.data if compound.data else {}
            logger.info(f"Updating compound '{compound.name}' with new data.")
            compound_data.update(new_data)
            compound.data = compound_data
            compound.save()
            logger.info(f"Successfully updated compound '{compound.name}'.")
        except Exception as e:
            logger.error(f"Error updating compound '{compound.name}': {str(e)}")

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
        if not input_file:
            logger.error("Input file is required.")
            return

        logger.info(f"Starting compound data import from {input_file}")
        self.import_compound_data(input_file)
        logger.info("Compound data import completed.")
