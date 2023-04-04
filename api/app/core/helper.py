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
