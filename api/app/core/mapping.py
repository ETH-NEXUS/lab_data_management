import csv
import re
from os.path import isfile
from typing import Union

from .helper import charToAlphaPos, posToAlphaChar


class Mapping:
    """Represents one mapping"""

    def __init__(
        self,
        from_pos: int,
        to_pos: int,
        amount: float = 0,
        status: str = None,
        map_type: bool = False,
        current_amount: float = 0,
        current_dmso: float = 0,
    ):
        self.__from = int(from_pos)
        self.__to = int(to_pos)
        self.__amount = float(amount)
        self.__status = status  # not str(status)
        self.__map_type = map_type  # if we need to map the well type (apply when we map a control plate)
        self.current_amount = current_amount
        self.current_dmso = current_dmso

    def __str__(self):
        return f"Mapping: {self.__from} -> {self.__to} ({self.__amount})"

    @property
    def from_pos(self) -> int:
        return self.__from

    @property
    def to_pos(self) -> int:
        return self.__to

    @property
    def amount(self) -> float:
        return self.__amount

    @property
    def status(self) -> str:
        return self.__status

    @property
    def map_type(self) -> bool:
        return self.__map_type


class MappingList:
    """
    Represents a list of mappings.
    This can be used to map one plate to another.
    """

    def __init__(self, source=None, target=None):
        self.__source = source
        self.__target = target
        self.__mappings = []
        self.current = 0

    def add(self, mapping: Mapping):
        self.__mappings.append(mapping)

    @property
    def source(self):
        return self.__source

    @property
    def target(self):
        return self.__target

    def __iter__(self):
        self.current = 0
        return self

    def __next__(self):
        if self.current < len(self.__mappings):
            ret = self.__mappings[self.current]
            self.current += 1
            return ret
        else:
            raise StopIteration

    def __getitem__(self, item):
        return self.__mappings[item]

    @classmethod
    def one_to_one(cls, n_pos: int, amount: float = 0, map_type: bool = False):
        """returns a one_to_one mapping for n_pos positions"""
        mappings = cls()
        for p in range(n_pos):
            print("adding mapping", p)
            mappings.add(Mapping(p, p, amount, map_type=map_type))
        return mappings

    @classmethod
    def from_csv(
        cls,
        csv_file: str,
        from_col: Union[str, int],
        to_col: Union[str, int],
        amount_col: Union[str, int],
        delimiter: str = ",",
        quotechar: str = '"',
    ) -> "MappingList":
        """
        Returns a MappingList generated from csv.
        The values are taken from the given colum names of indexes.
        """
        if isfile(csv_file):
            mappings = cls()
            with open(csv_file, "r", newline="") as cf:
                if isinstance(from_col, str):
                    reader = csv.DictReader(
                        cf, delimiter=delimiter, quotechar=quotechar
                    )
                else:
                    reader = csv.reader(cf, delimiter=delimiter, quotechar=quotechar)
                for row in reader:
                    mappings.add(Mapping(row[from_col], row[to_col], row[amount_col]))
            return mappings
        else:
            raise FileNotFoundError(f"Cannot find csv file {csv_file}.")


class PositionMappingError:
    def __init__(self, position, message="Cannot convert position to row, col: {}"):
        self.message = message.format(position)
        super().__init__(self.message)


class PositionMapper:
    """
    Maps a string notation position to row, col and the other way round.
    i.e. A2 -> row: 0, col: 1
    """

    @staticmethod
    def map(position: str) -> tuple[int, int]:
        match = re.match(r"(?P<alpha>[A-Z]+)(?P<num>[0-9]+)", position, re.IGNORECASE)
        if match:
            row = charToAlphaPos(match.group("alpha"))
            col = int(match.group("num"))

        else:
            raise PositionMappingError(position)
        return row, col

    @staticmethod
    def unmap(row: int, col: int) -> str:
        return f"{posToAlphaChar(row)}{col}"


# @staticmethod
#     def convert_position_to_index(position: str,
#                                   number_of_columns: int) -> int:
#         match = re.match(r'([a-zA-Z]+)(\d+)', position)
#         letters = match.group(1)
#         col = int(match.group(2))
#         row = 0
#         for index, char in enumerate(letters[::-1]):
#             row += (ord(char.upper()) - ord('A') + 1) * (26 ** index)
#         index = (row - 1) * number_of_columns + (col - 1)
#         return index
