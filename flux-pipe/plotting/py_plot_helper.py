"""_summary_

Returns:
    _type_: _description_
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


## Helper functions to parse the output file
def convert_value(value):
    """ Convert a string to an int or float if possible, otherwise return the string.

    Parameters
    ----------
    value : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value.strip()

def parse_line(line):
    """ Parse a line of the output file into a dictionary.

    Parameters
    ----------
    line : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """

    key_values = line.strip().split(',')
    # print(key_values)
    parsed_dict = {}
    for key_value in key_values:
        if ':' in key_value:
            key, value = key_value.split(':')
            parsed_dict[key.strip()] = convert_value(value)
    return parsed_dict

def load_data(_file_path):
    """Load the data from the file into a pandas DataFrame.

    Parameters
    ----------
    _file_path : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """

    data = []
    print("\n", _file_path, '\n')
    with open(_file_path, 'r', encoding="utf-8") as file:
        for line in file.readlines():
            data.append(parse_line(line))
    return pd.DataFrame(data)


def get_ax(ax=None):
    """ Get the fig and ax, either a new one or the one passed in.

    Parameters
    ----------
    ax : _type_, optional
        _description_, by default None

    Returns
    -------
    _type_
        _description_
    """
    if ax is not None:
        fig, ax0 = ax.get_figure(), ax
    else:
        fig, ax0 = plt.subplots(1)
    return fig, ax0

def add_parent_dir(): #depricated and wrong
    """_summary_
    """
    # Get the parent directory path
    this_cwd = os.getcwd()
    if "plotting" in this_cwd:
        pipdir = os.path.abspath(os.path.join(this_cwd, ".."))
    elif "pipe" in this_cwd:
        pipdir = this_cwd
    else:
        pipdir = os.path.abspath(os.path.join(this_cwd, "flux-pipe"))

    # Add the parent directory to the module search path
    sys.path.append(pipdir)

def add_pipedir(do_print=False):
    """_summary_

    Parameters
    ----------
    do_print : bool, optional
        _description_, by default False
    """

    # get the pipedir environment variable
    pipedir = os.environ.get("PIPEPATH")
    #add the directory to the search path
    sys.path.append(pipedir)
    if do_print:
        print(f"\nAdded {pipedir} to the module search path!\n")

# Remove outliers from the dataset
def remove_outliers(data, ph0_temp, th0_temp, threshold=3):
    """_summary_

    Parameters
    ----------
    data : _type_
        _description_
    ph0_temp : _type_
        _description_
    th0_temp : _type_
        _description_
    threshold : int, optional
        _description_, by default 3

    Returns
    -------
    _type_
        _description_
    """

    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    filtered_data, good_points = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data)
                                    if mean - threshold * std < x < mean + threshold * std]))
    ph0c, th0c = ph0_temp[good_points], th0_temp[good_points]
    bad_points = [i for i in range(len(data)) if i not in good_points]
    ph0_b, th0_b = ph0_temp[bad_points], th0_temp[bad_points]
    outlier_data = data[bad_points]
    mean = np.mean(filtered_data)
    std = np.std(filtered_data)
    return np.asarray(filtered_data), mean, std, ph0c, th0c, outlier_data, ph0_b, th0_b


def scale_data(vel0_clean, vel1_clean, outlier_V0, outlier_V1, scale=15**2, power=1):
    """ Scale the data between 0 and 1, then raise to a power, then scale by a factor

    Parameters
    ----------
    vel0_clean : _type_
        _description_
    vel1_clean : _type_
        _description_
    outlier_V0 : _type_
        _description_
    outlier_V1 : _type_
        _description_
    scale : _type_, optional
        _description_, by default 15**2
    power : int, optional
        _description_, by default 1

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    to
        _description_
    """

    vel0_max = np.nanmax(vel0_clean)
    vel1_max = np.nanmax(vel1_clean)
    vel0_min = np.nanmin(vel0_clean)
    vel1_min = np.nanmin(vel1_clean)

    v0 = scale * ((np.abs(vel0_clean) - vel0_min) / (vel0_max-vel0_min))**power
    v1 = scale * ((np.abs(vel1_clean) - vel1_min) / (vel1_max-vel1_min))**power

    outlier_V0_scaled = scale * ((np.abs(outlier_V0) - vel0_min) / (vel0_max-vel0_min))**power
    outlier_V1_scaled = scale * ((np.abs(outlier_V1) - vel1_min) / (vel1_max-vel1_min))**power

    return v0, v1, outlier_V0_scaled, outlier_V1_scaled


# add_pipedir()
