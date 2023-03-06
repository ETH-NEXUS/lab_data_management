from django.utils.translation import gettext as _
from os import environ
from rest_framework.test import APIClient
import pandas as pd
import numpy as np
from typing import Callable
from core.models import Plate, Well
from core.helper import posToAlphaChar


api = APIClient()


def login(username: str = None, password: str = None):
    """ Logs in to the api. """
    from getpass import getpass
    if not username:
        username = input('Username:')
    if not password:
        password = getpass(prompt='Password:')
    if api.login(username=username, password=password):
        print("Logged in.")
    else:
        print("Login failed!")


def version():
    """ Returns the git version of LDM. """
    print(f"LDM Version: {environ.get('GIT_VERSION', 'N/A')}")


def df(url: str):
    """ Returns a pandas DataFrame from the api query. """
    resp = api.get(url)
    if resp and resp.status_code == 200:
        if 'results' in resp.data:
            return pd.DataFrame.from_dict(resp.data.get('results'))
        else:
            return pd.DataFrame.from_dict(resp.data)


def col2df(df: pd.DataFrame, col: str):
    """
    Return a Pandas dataframe based on a column that contains a list of dicts.
    """
    rows = []
    for _, row in df[col].items():
        for item in row:
            rows.append(item)
    df = pd.DataFrame(rows)
    return df


def find(lst: list, _filter: Callable):
    """
    returns the first element that matches the filter
    """
    for elem in lst:
        if _filter(elem):
            return elem
    return None


def plate(barcode: str):
    """ 
    Returns a pandas DataFrame with the wells of the 
    plate with the given barcode
    """
    return df(f"/api/plates/?barcode={barcode}")


def measurement(barcode: str, abbrev: str, matrix: bool = False):
    """ 
    Returns a pandas DataFrame with the measurement values given by
    abbrev. The DataFrame can be used to display a heatmap of the 
    plate using:
    ```
    library('bioassays')
    matrix = matrix96(df, 'value', rm='TRUE')
    heatplate(matrix, "Plate", size=15)
    ```
    """
    plate = Plate.objects.get(barcode=barcode)
    if matrix:
        df = pd.DataFrame(np.nan, index=[posToAlphaChar(row) for row in range(1, plate.dimension.rows + 1)], columns=range(1, plate.dimension.cols + 1))
    else:
        df = pd.DataFrame(np.nan, index=range(plate.num_wells), columns=('row', 'col', 'position', 'value'))

    for position in range(plate.num_wells):
        row, col = plate.dimension.row_col(position)
        row = posToAlphaChar(row)
        try:
            well = Well.objects.get(plate=plate, position=position)
            if matrix:
                df.loc[row, col] = well.measurement(abbrev) or np.nan
            else:
                df.iloc[position] = [row, col, position, well.measurement(abbrev) or np.nan]
        except Well.DoesNotExist:
            if matrix:
                df.loc[row, col] = np.nan
            else:
                df.iloc[position] = [row, col, position, np.nan]
    return df
