"""
Helper Functions for Plotting Data from Fluxon Simulations
==========================================================

Description
-----------
This module contains a collection of helper functions designed to assist in parsing and
plotting data generated from fluxon simulations.

Functions
---------
    convert_value(value: str) -> Union[int, float, str]:
        Convert a string to an int or float if possible, otherwise return the string.

    parse_line(line: str) -> Dict[str, Union[int, float, str]]:
        Parse a line of the output file into a dictionary.

    load_data(the_path: str) -> pd.DataFrame:
        Load the data from the file into a pandas DataFrame.

    get_ax(ax: Optional[matplotlib.axis]) -> Tuple[matplotlib.figure, matplotlib.axis]:
        Get or create a pyplot figure and axis pair.

    add_fluxon_dirs_to_path(do_print: bool = False) -> None:
        Add the fluxon directories to the system path.

    list_directories(path: str) -> List[str]:
        List the directories in the given path.

Modules
-------
    os :
        For directory and file operations.

    sys :
        For system-specific parameters and functions.

    matplotlib.pyplot :
        For plotting.

    numpy :
        For numerical operations.

    pandas :
        For data manipulation and analysis.

Example
-------
    ```python
    # Example usage of convert_value function
    result = convert_value("42")
    ```

Raises
------
    FileNotFoundError:
        If the specified file in `load_data` cannot be found.

    ValueError:
        If the data file cannot be loaded or parsed correctly.

See Also
--------
    Other related scripts and modules for fluxon simulations.

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
    value : string
        the value to convert to a number

    Returns
    -------
    int, float, string
        the input value, typecast appropriately
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
    line : str
        one line of the output file, with key:value pairs separated by commas

    Returns
    -------
    dict
        a dictionary of the key:value pairs
    """

    key_values = line.strip().split(',')
    # print(key_values)
    parsed_dict = {}
    for key_value in key_values:
        if ':' in key_value:
            key, value = key_value.split(':')
            parsed_dict[key.strip()] = convert_value(value)
    return parsed_dict


def load_data(the_path):
    """Load the data from the file into a pandas DataFrame.

    Parameters
    ----------
    the_path : str
        the path to the file to load

    Returns
    -------
    dataframe
        the data from the file given by the_path
    """

    data = []
    print("\n", the_path, '\n')
    with open(the_path, 'r', encoding="utf-8") as file:
        for line in file.readlines():
            data.append(parse_line(line))
    return pd.DataFrame(data)


def get_ax(ax=None):
    """ Get the fig and ax. If None, create a new fig and ax.
    Otherwise, return the given ax.

    Parameters
    ----------
    ax : pyplot axis or None, optional
        Either an axis or None, by default None

    Returns
    -------
    figure, axis
        a pyplot figure and axis pair
    """
    if ax is not None:
        fig, ax0 = ax.get_figure(), ax
    else:
        fig, ax0 = plt.subplots(1)
    return fig, ax0


def add_fluxon_dirs_to_path(do_print=False):
    """ Add the fluxon directories to the system path.

    Parameters
    ----------
    do_print : bool, optional
        print the paths added, by default False
    """

    # Get the current directory path
    this_cwd = os.getcwd()

    # Get the list of directories in the current directory
    dirs = list_directories(this_cwd)
    dirlist = [os.path.join(this_cwd, x) for x in dirs if "fluxon" in x]

    # Add the pipe and plotting directories to the path
    for thepath in dirlist:
        if "mhd" in thepath:
            dirlist.append(os.path.join(thepath, "flux-pipe"))
            dirlist.append(os.path.join(thepath, "flux-pipe", "plotting"))
            break

    # Get the pipedir environment variable and add it to the path
    pipedir = os.environ.get("PIPEPATH")
    if pipedir is not None:
        dirlist.append(pipedir)

    # Add the parent directory to the module search path
    for path in dirlist:
        sys.path.append(path)
        if do_print:
            print(path)


def list_directories(path):
    """ List the directories in the given path.

    Parameters
    ----------
    path : str
        the directory to list the subdirectories of

    Returns
    -------
    list
        a list of the subdirectories of the given path
    """

    dirs = []
    with os.scandir(path) as entries:
        for entry in entries:
            if entry.is_dir():
                dirs.append(entry.name)
    return dirs


add_fluxon_dirs_to_path()
