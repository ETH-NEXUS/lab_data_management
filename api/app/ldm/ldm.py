from os import environ
from rest_framework.test import APIClient
import pandas as pd
import numpy as np
from typing import Callable
from core.models import Plate, Well, PlateDetail
from core.helper import posToAlphaChar

# from scipy.stats import median_abs_deviation as mad

api = APIClient()


def login(username: str = None, password: str = None):
    """Logs in to the api."""
    from getpass import getpass

    if not username:
        username = input("Username:")
    if not password:
        password = getpass(prompt="Password:")
    if api.login(username=username, password=password):
        print("Logged in.")
    else:
        print("Login failed!")


def version():
    """Returns the git version of LDM."""
    print(f"LDM Version: {environ.get('GIT_VERSION', 'N/A')}")


def df(url: str):
    """Returns a pandas DataFrame from the api query."""
    resp = api.get(url)
    if resp and resp.status_code == 200:
        if "results" in resp.data:
            return pd.DataFrame.from_dict(resp.data.get("results"))
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
        df = pd.DataFrame(np.nan, index=[posToAlphaChar(row) for row in range(1, plate.dimension.rows + 1)], columns=range(1, plate.dimension.cols + 1), )
    else:
        df = pd.DataFrame(np.nan, index=range(plate.num_wells), columns=("row", "col", "position", "value"), )

    for position in range(plate.num_wells):
        row, col = plate.dimension.row_col(position)
        row = posToAlphaChar(row)
        try:
            well = Well.objects.get(plate=plate, position=position)
            if matrix:
                df.loc[row, col] = well.measurement(abbrev) or np.nan
            else:
                df.iloc[position] = [row, col, position, well.measurement(abbrev) or np.nan, ]
        except Well.DoesNotExist:
            if matrix:
                df.loc[row, col] = np.nan
            else:
                df.iloc[position] = [row, col, position, np.nan]
    return df


def get_experiment_measurements(experiment_name: str):
    """
    Returns a pd DataFrame of measurements for a given experiment.
    """
    experiment_plates = Plate.objects.filter(experiment__name=experiment_name)

    rows = []
    for pl in experiment_plates:
        plate_dimension = pl.dimension
        print(f"Processing plate {pl.barcode}")

        wells = pl.wells.all()  # Well.objects.filter(plate=pl)
        for well in wells:
            row, col = plate_dimension.row_col(well.position)

            measurements = well.measurements.all()
            well_rows = [{
                "WellCoordinate": well.hr_position, "Value": measurement.value, "Plate": pl.barcode, "plateRow": row, "plateColumn": col, "Control": well.type.name,
                "measurement": measurement.label,
            } for measurement in measurements]
            rows.extend(well_rows)

    rows = sorted(rows, key=lambda k: k["measurement"])

    print(f"\nFound {len(rows)} measurements")
    return pd.DataFrame(rows)


def normalize_values(raw_data, log_value=False, label=None, pos_neg_only=False):
    raw_data = raw_data.copy()

    if pos_neg_only:
        raw_data = raw_data[raw_data['Control'].isin(['P', 'N'])]
    if label is not None:
        raw_data = raw_data[raw_data['measurement'] == label]
    if log_value:
        raw_data['Value'] = np.log10(raw_data['Value'])

    plate_median = raw_data.groupby(['Plate', 'Control'])['Value'].median().reset_index()
    plate_median = plate_median.pivot_table(index='Plate', columns='Control', values='Value').reset_index()

    column_names_map = {'C': 'med.c', 'N': 'med.n', 'P': 'med.p'}
    plate_median.columns = ['Plate'] + [column_names_map[c] for c in plate_median.columns[1:]]

    result = pd.merge(raw_data, plate_median, on=['Plate'])
    result['norm'] = (result['Value'] - result['med.n']) / (result['med.p'] - result['med.n'])

    return result


def mad(series, constant=1.4826):
    # this is how it is done in R, it uses this scaling factor by default (the median_abs_deviation from spy.stats gives different results)
    # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/mad
    return constant * (series - series.median()).abs().median()


def calculate_z_prime(raw_data, log_value=False):
    data = raw_data.copy()

    if log_value:
        data['Value'] = np.log10(data['Value'])

    plate_med = data.groupby(['Plate', 'Control'])['Value'].median().unstack().reset_index()
    plate_mad = data.groupby(['Plate', 'Control'])['Value'].apply(mad).unstack().reset_index()

    z_prime = 1 - (3 * (plate_mad['P'] + plate_mad['N']) / abs(plate_med['P'] - plate_med['N']))

    result = pd.concat([plate_med['Plate'], z_prime], axis=1)
    result.columns = ['Plate', 'z_prime']

    return result
