
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
