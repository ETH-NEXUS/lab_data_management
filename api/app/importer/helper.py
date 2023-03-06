
def sameSchema(dict1, dict2, same=True):
    """
    Compares the keys of two dicts and returns true
    if they are identical false otherwise
    """

    for key in dict1.keys():
        if key in dict2:
            if isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
                same &= sameSchema(dict1[key], dict2[key])
            else:
                same &= True
        else:
            return False
    return same


plate_dimension_lookup = {
    96: (8, 12),
    384: (16, 24),
    1536: (32, 48)
}


def row_col_from_wells(wells: int):
    return plate_dimension_lookup[wells]


def closest(val: int, values: tuple[int]):
    """
    Returns the closest number in a list of numbers.
    If the value is exactly between the smaller value is taken.
    """
    return values[min(range(len(values)), key=lambda i: abs(values[i]-val))]


def normalize_row(row: int):
    return closest(row, (8, 16, 32))


def normalize_col(col: int):
    return closest(col, (12, 24, 48))
