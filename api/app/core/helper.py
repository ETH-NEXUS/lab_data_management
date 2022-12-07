from string import ascii_lowercase


LETTERS = {letter: str(index)
           for index, letter in enumerate(ascii_lowercase, start=1)}


def charToAlphaPos(chr: str):
    """
    Maps a character to a number
    A,a -> 1
    B,b -> 2
    ...
    Z,z -> 26
    """
    if len(chr) != 1:
        raise Exception("Can only convert one char")
    return int(LETTERS[chr.lower()])


def posToAlphaChar(pos: int):
    """
    Maps a number to a character
    1 -> A
    2 -> B
    ...
    26 -> Z
    """
    return ascii_lowercase[pos - 1].upper()
