from .helper import charToAlphaPos, posToAlphaChar
import re


class Mapping():
    """Represents one mapping"""

    def __init__(self, from_pos: int, to_pos: int, amount: float = 0):
        self.__from = from_pos
        self.__to = to_pos
        self.__amount = amount

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

    @classmethod
    def one_to_one(cls, n_pos: int, amount: float = 0):
        """returns a one_to_one mapping for n_pos positions"""
        mappings = cls()
        for p in range(n_pos):
            mappings.add(Mapping(p, p, amount))
        return mappings


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
