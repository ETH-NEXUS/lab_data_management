"""
Lists of well-to-well transfers, used to map one plate onto another (Plate.map).
A list is built one to one for a plate copy, from a csv file, or by the importer.
"""

import csv
from os.path import isfile
from typing import Union


class Mapping:
    """Represents one mapping"""

    def __init__(
        self,
        from_pos: int,
        to_pos: int,
        amount: float = 0,
        status: str | None = None,
        map_type: bool = False,
        # A mapping that carries no instrument reading leaves these at None,
        # which means "not reported". Zero would mean "the well is empty".
        current_amount: float | None = None,
        current_dmso: float | None = None,
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
    def status(self) -> str | None:
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
