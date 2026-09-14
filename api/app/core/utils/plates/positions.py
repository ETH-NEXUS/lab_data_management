"""
Well positions on a plate: the letters of a row and the "A3" notation.

Example: "A3" <-> row 1, column 3, and row letters A -> 1, Z -> 26, AA -> 27.
"""

from string import ascii_uppercase
import re


def charToAlphaPos(letters: str):
    """
    Maps a character sequence to a number
    A,a -> 1
    B,b -> 2
    ...
    Z,z -> 26
    AA,aa -> 27
    AB,ab -> 28
    ...
    AZ,az -> 52
    """
    if not re.match(r"[A-z]+", letters):
        raise ValueError("Only letters are allowed!")
    pos = 0
    for index, char in enumerate(letters[::-1]):
        pos += (ord(char.upper()) - ord("A") + 1) * (26**index)
    return pos


def posToAlphaChar(pos: int):
    """
    Maps a number to a character
    1 -> A
    2 -> B
    ...
    26 -> Z
    ...
    27 -> AA
    28 -> AB
    ...
    52 -> AZ
    """
    try:
        letter = ""
        while pos > 0:
            pos, remainder = divmod(pos - 1, 26)
            letter = ascii_uppercase[remainder] + letter
        return letter
    except IndexError:
        return "?"


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
