from os.path import isfile
import re
import csv

from .helper import charToAlphaPos, posToAlphaChar


class Mapping():
    """Represents one mapping"""

    def __init__(self, from_pos: int, to_pos: int, amount: float = 0):
        self.__from = int(from_pos)
        self.__to = int(to_pos)
        self.__amount = float(amount)

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


class MappingList():
    """
    Represents a list of mappings.
    This can be used to map one plate to another.
    """

    def __init__(self):
        self.__mappings = []
        self.current = 0

    def add(self, mapping: Mapping):
        self.__mappings.append(mapping)

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
    def one_to_one(cls, n_pos: int, amount: float = 0):
        """returns a one_to_one mapping for n_pos positions"""
        mappings = cls()
        for p in range(n_pos):
            mappings.add(Mapping(p, p, amount))
        return mappings

    @classmethod
    def from_csv(cls, csv_file: str, from_col: str | int, to_col: str | int, amount_col: str | int, delimiter: str = ',', quotechar: str = '"') -> 'MappingList':
        """
        Returns a MappingList generated from csv. 
        The values are taken from the given colum names of indexes.
        """
        if isfile(csv_file):
            mappings = cls()
            with open(csv_file, 'r', newline='') as cf:
                if isinstance(from_col, str):
                    reader = csv.DictReader(cf, delimiter=delimiter, quotechar=quotechar)
                else:
                    reader = csv.reader(cf, delimiter=delimiter, quotechar=quotechar)
                for row in reader:
                    mappings.add(Mapping(
                        row[from_col],
                        row[to_col],
                        row[amount_col]
                    ))
            return mappings
        else:
            raise FileNotFoundError(f"Cannot find csv file {csv_file}.")


class PositionMappingError():
    def __init__(self, position, message="Cannot convert position to row, col: {}"):
        self.message = message.format(position)
        super().__init__(self.message)


class PositionMapper():
    """ 
    Maps a string notation position to row, col and the other way round. 
    i.e. A2 -> row: 0, col: 1
    """
    @staticmethod
    def map(position: str) -> tuple[int, int]:
        match = re.match(
            r'(?P<alpha>[A-Z])(?P<num>[0-9]+)', position, re.IGNORECASE)
        if match:
            row = charToAlphaPos(match.group('alpha'))
            col = int(match.group('num'))
        else:
            raise PositionMappingError(position)
        return row, col

    @staticmethod
    def unmap(row: int, col: int) -> str:
        return (f"{posToAlphaChar(row)}{col}")
