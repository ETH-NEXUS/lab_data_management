from typing import List
from helpers.logger import logger


class CompoundInfoReader:
    def __init__(self, column_names_file_path="compound_info_columns.txt"):
        self.column_names_file_path = column_names_file_path

    def read(self) -> List:
        with open(self.column_names_file_path, "r") as file:
            logger.info(f"Reading column names from {self.column_names_file_path}")
            column_names = [name.strip().lstrip() for name in file.readlines()]
            logger.debug(f"Column names: {column_names}")
            return column_names
