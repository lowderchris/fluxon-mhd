"""
pipe_helper: Comprehensive Library for Flux-Pipe Algorithm and Fluxon Simulations
===============================================================================

This library provides a collection of utility functions to assist with the
Flux-Pipe algorithm and Fluxon simulations. It offers functionalities for managing directories,
handling FITS files, manipulating magnetogram data, parsing and plotting data generated from fluxon simulations.

Modules:
--------
- os, os.path, sys, pathlib.PosixPath, Path
- numpy as np
- matplotlib.pyplot as plt
- astropy.io.fits, astropy.nddata.block_reduce
- sunpy, sunpy.coordinates, sunpy.io, sunpy.net.Fido, attrs as a
- pandas

Functions:
----------
### General Utilities
- `configurations`: Reads and sanitizes configuration settings from a specified config file.
- `convert_value`: Converts a string to an int or float if possible.
- `calculate_directories`: Helper function to calculate directories.
- `add_dir_to_path`: Adds a directory and all subdirectories to the PATH environment variable.
- `add_top_level_dirs_to_path`: Adds the top-level directories under a root directory to the PATH.
- `add_paths`: Adds various paths to the system path.
- `find_file_with_string`: Searches a directory for a file containing a given string.
- `shorten_path`: Removes the DATAPATH environment variable from a string.

### Magnetogram Utilities
- `make_mag_dir`: Creates a directory for magnetogram data.
- `get_magnetogram_file`: Grabs HMI data.
- `reduce_mag_file`: Reduces the size of a magnetogram FITS file by a given factor.
- `reduce_fits_image`: Reduces the size of a FITS image.
- `plot_raw_magnetogram`: Plots the magnetogram.
- `load_fits_magnetogram`: Loads a magnetogram from a FITS file.
- `write_magnetogram_params`: Writes the magnetic_target.params file for a given CR and reduction amount.
- `load_magnetogram_params`: Reads the magnetic_target.params file and returns the parameters.
- `read_fits_data`: Reads FITS data and fixes/ignores any non-standard FITS keywords.
- `get_fixed_coords`: Corrects input coordinates.

### Fluxon Simulation Utilities
- `parse_line`: Parse a line of the output file into a dictionary.
- `load_data`: Load the data from the file into a pandas DataFrame.
- `get_ax`: Get or create a pyplot figure and axis pair.
- `add_fluxon_dirs_to_path`: Add the fluxon directories to the system path.
- `list_directories`: List the directories in the given path.
- `path_add`: Add directories to the system path.

Usage Example:
--------------
```python
# Example usage of convert_value function
import pipe_helper as ph
result = ph.convert_value("42")

Author:
-------
    Gilly <gilly@swri.org> (and others!)

Dependencies:
-------------
    os, os.path, sys, pathlib.PosixPath, Path, numpy as np,
    matplotlib.pyplot as plt, astropy.io.fits, astropy.nddata.block_reduce,
    sunpy, sunpy.coordinates, sunpy.io, sunpy.net.Fido, attrs as a, pandas
"""


# Import libraries
import os
import os.path
import sys
import ast
from pathlib import PosixPath, Path

import pandas as pd
# from pipe_helper import convert_value
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.nddata import block_reduce

import sunpy
import sunpy.coordinates
import sunpy.io
from sunpy.net import Fido, attrs as a
import configparser



def configurations(config_name=None, config_filename="config.ini", debug=False):
    """
    Reads and sanitizes configuration settings from a specified config file.

    Args:
        config_name (str, optional): The specific configuration section to read.
                                     If not provided, defaults to the DEFAULT section.
        config_filename (str, optional): The filename of the configuration file.
                                         Defaults to 'config.ini'.
        debug (bool, optional): Whether to print debug information. Defaults to False.

    Returns:
        dict: Configuration settings as key-value pairs.
    """
    config_obj = configparser.ConfigParser()
    config_path = f"fluxon-mhd/flux-pipe/config/{config_filename}"

    # Search for the configuration file in the current directory and subdirectories
    if not os.path.exists(config_path):
        found = False
        for root, dirs, files in os.walk(os.getcwd()):
            if config_filename in files:
                config_path = os.path.join(root, config_filename)
                found = True
                break
        if not found:
            raise FileNotFoundError("Configuration file not found.")

    # Clean the file content: remove comments and trailing whitespaces
    with open(config_path, "r") as f:
        lines = f.readlines()
    clean_lines = [line.split("#")[0].rstrip() for line in lines]
    clean_content = "\n".join(clean_lines)

    # Parse the clean configuration string
    config_obj.read_string(clean_content)

    # Fallback to section defined in the DEFAULT section if no specific section is provided
    if config_name is None:
        config_name = config_obj["DEFAULT"]['config_name']

    # Extract and further process configuration settings
    the_config = dict(config_obj[config_name])
    the_config['abs_rc_path'] = os.path.expanduser(the_config['rc_path'])
    the_config["run_script"] = os.path.join(the_config['fl_prefix'], the_config["run_script"])
    the_config["rotations"] = ast.literal_eval(the_config["rotations"])
    the_config["fluxon_count"] = ast.literal_eval(the_config["fluxon_count"])
    the_config["n_jobs"] = str(len(the_config["rotations"]) * len(the_config["fluxon_count"]))

    # Calculate directories
    calculate_directories(the_config)

    the_config['magfile'] = f"{the_config['mag_dir']}/CR{{}}_r{the_config['mag_reduce']}_hmi.fits"

    for key, value in the_config.items():
        the_config[key] = convert_value(value)

    if debug:
        print("\nConfiguration file values:\n--------------------------------")
        for key, value in sorted(the_config.items()):
            print(f"{key}: \t{value}")
        print("--------------------------------\n\n")

    return the_config


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

    if type(value) is list:
        new_list = []
        for item in value:
            new_list.append(convert_value(item))
        return new_list

    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value.strip()

# Helper function to calculate directories
def calculate_directories(the_config):
    basedir = the_config['base_dir'].strip()
    batch_name = the_config['batch_name'].strip()
    dat_dir = the_config.get('data_dir', None)  # Assuming you have this in your config

    fluxdir = os.path.join(basedir, "fluxon-mhd")
    pipedir = os.path.join(fluxdir, "flux-pipe")
    pdldir = os.path.join(fluxdir, "pdl", "PDL")

    # Use the provided data_dir if defined, otherwise calculate it
    datdir = dat_dir if dat_dir else os.path.join(basedir, "fluxon-data")

    magdir = os.path.join(datdir, "magnetograms")
    batchdir = os.path.join(datdir, "batches", batch_name)
    logfile = os.path.join(batchdir, "pipe_log.txt")

    # Update the original config dictionary
    the_config['pipe_dir']  = pipedir
    the_config['pdl_dir']   = pdldir
    the_config['datdir']    = datdir
    the_config['mag_dir']   = magdir
    the_config['batch_dir'] = batchdir
    the_config['logfile']   = logfile


configs = configurations()
dat_dir = configs["data_dir"]
default_email = configs["jsoc_email"]

# PATH MANAGEMENT

def add_dir_to_path(root_dir=None):
    """Adds a directory and all subdirectories to the PATH environment variable.

    Parameters
    ----------
    root_dir : str, optional
        Root directory path
    """

    if root_dir is None:
        root_dir = os.environ("FL_PREFIX", None) or "fluxon-mhd/"

    # Get the current PATH
    current_path = os.environ.get('PATH', '')

    # Initialize a set with the current PATH elements to avoid duplicates
    path_set = set(current_path.split(os.pathsep))

    # Walk through the directory tree
    for dirpath, _, _ in os.walk(root_dir):
        # Add each directory to the set
        path_set.add(dirpath)

    # Convert the set back to a string
    new_path = os.pathsep.join(path_set)

    # Update the PATH
    os.environ['PATH'] = new_path


def add_top_level_dirs_to_path(root_dir):
    """Adds the top-level directories under a root directory to the PATH environment variable.

    Parameters
    ----------
    root_dir : str
        Root directory path
    """

    # Get the current PATH
    current_path = os.environ.get('PATH', '')

    # Initialize a set with the current PATH elements to avoid duplicates
    path_set = set(current_path.split(os.pathsep))

    # List the top-level directories under the root directory
    top_level_dirs = [os.path.join(root_dir, d) for d in os.listdir(root_dir) \
                      if os.path.isdir(os.path.join(root_dir, d))]
    top_level_dirs.append(f"{root_dir}:")

    # Add each top-level directory to the set
    path_set.update(top_level_dirs)

    # Convert the set back to a string
    new_path = os.pathsep.join(path_set)

    # Update the PATH
    os.environ['PATH'] = new_path

    # news = os.environ.get('PATH', '')
    # news_list = news.split(os.pathsep)
    # news_list.sort()
    # a=[print(x) for x in news_list]
    return new_path


def set_paths(do_plot=False):
    """
    Checks if the environment variable FL_PREFIX is set and optionally prints the paths.

    Parameters:
        do_plot (bool): Whether to print the paths or not.

    Returns:
        None
    """

    fl_prefix = os.environ.get('FL_PREFIX')

    if fl_prefix:
        envpath = os.path.join(fl_prefix, 'python_paths.py')

        # Check if the file exists and is readable
        if os.path.exists(envpath) and os.access(envpath, os.R_OK):
            exec(open(envpath).read())
        else:
            print(f"File does not exist or is not readable: {envpath}")
    else:
        print("Environment variable FL_PREFIX is not set.")

    # Optionally print the lists of directories
    if do_plot:
        print("\n\nsys.path has:")
        for path in sys.path:
            print(f" {path}")
        print("--------------------------\n")

        # Add any additional paths you want to print here
        # For example, if you have a list similar to Python's sys.path
        # print("\nYour_Path has:")
        # for path in Your_Path:
        #     print(f" {path}")
        # print("--------------------------\n")


def add_paths(flux_pipe_dir):
    """Adds various paths to the system path.

    Parameters
    ----------
    flux_pipe_dir : str
        FLUXpipe directory path

    Returns
    -------
    str
        PDL script directory path
    """

    # Path to the PDL script
    pdl_script_path = flux_pipe_dir + "magnetogram2wind.pdl"
    os.chdir(flux_pipe_dir)
    # Get the plotscript directory path
    plot_dir = os.path.abspath(os.path.join(flux_pipe_dir, "plotting"))
    sys.path.append(plot_dir)
    return pdl_script_path


def find_file_with_string(directory, search_string):
    """Searches a directory for a file containing a given string.

    Parameters
    ----------
    directory : str
        Search directory path
    search_string : str
        Search string

    Returns
    -------
    str
        Search result file path

    """
    for file_name in os.listdir(directory):
        if search_string in file_name:
            return os.path.join(directory, file_name)
    return None


def shorten_path(string):
    """Removes the DATAPATH environment variable from a string.
    This makes it much more readable when printing paths.

    Parameters
    ----------
    string : str
        String to shorten

    Returns
    -------
    str
        Shortened string
    """
    datapath = os.getenv("DATAPATH")
    if datapath:
        return string.replace(datapath, "$DATAPATH")
    else:
        return string


# MAGNETOGRAM MANAGEMENT


def make_mag_dir(datdir):
    """Creates a directory for magnetogram data.

    Parameters
    ----------
    datdir : str
        Data directory path

    Returns
    -------
    str
        Magnetogram directory path
    """

    mag_dir = os.path.join(datdir, "magnetograms")
    if not os.path.exists(mag_dir):
        os.makedirs(mag_dir)
    return mag_dir


def get_magnetogram_file(cr=None, date=None, datdir=None, email=None,
                         force_download=False, reduce = False):
    """
    Grabs HMI data.

    Parameters
    ----------
    cr : int
        Carrington rotation number.
    date : str
        Date in YYYY-MM-DD format.
    data_dir : str, optional
        Directory where data will be stored. If not specified, default directories will be used.

    Returns
    -------
    big_path : str
        Path to the full resolution magnetogram file.
    small_path : str
        Path to the reduced resolution magnetogram file.
    """

    # Set the download account
    try:
        jsoc_email = email or os.environ["JSOC_EMAIL"]
    except KeyError:
        jsoc_email = default_email

    # Set the Carrington rotation number
    if cr is not None:
        CR = cr
        date = sunpy.coordinates.sun.carrington_rotation_time(CR)
    elif date is not None:
        CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))
    else:
        raise ValueError("Must specify either cr or date!")
    mag_dir = make_mag_dir(datdir)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"(py) Getting Magnetogram for CR{CR}, from {date}...")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    hmi_object = Path(mag_dir)
    file_list = list(hmi_object.iterdir())
    print("\tSearching for file...")
    found_file = False
    file = None
    for file in file_list:
        if str(CR)+"_r1_" in str(file):
            print(f"\t\tFound '{os.path.basename(file)}' in '{shorten_path(mag_dir)}'")
            found_file = True
            break

    if found_file:
        if force_download:
            print("\tForcing redownload...")
        else:
            small_path = reduce_mag_file(file, reduce, force=force_download)
            return file, small_path

    # c = drms.Client()
    # Generate a search
    crot = a.jsoc.PrimeKey('CAR_ROT', str(CR))
    res = Fido.search(a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), crot,
                    a.jsoc.Notify(jsoc_email))

    # Once the query is made and trimmed down...
    big_path = os.path.join(mag_dir, f"CR{CR}_r1_hmi.fits")
    # hmi_path = hmidat+'/{file}.fits'

    print("\tDownloading HMI from JSOC...")
    out = Fido.fetch(res, path=mag_dir)
    hmi_path_out = out[0]
    os.rename(hmi_path_out, big_path)
    print(f"\n\tSaved to {big_path}\n")

    small_path = reduce_mag_file(big_path, reduce, force=force_download)
    return big_path, small_path


def reduce_mag_file(mag_file, reduction=3, force=False):
    """Reduces the size of a magnetogram FITS file by a given factor.

    Parameters
    ----------
    mag_file : str
        Path to input magnetogram file
    reduction : int, optional
        Reduction factor. Defaults to 3.
    force : bool, optional
        Overwrite toggle. Defaults to False.

    Returns
    -------
    small_file: str
        Path to the reduced resolution magnetogram file.
    """

    small_file = PosixPath(str(mag_file).replace("_r1_", f"_r{reduction}_"))
    # reduce the FITS image
    print(f"\tReducing image size by a factor of {reduction}...", end="")
    if not os.path.exists(small_file) or force:
        small_file = reduce_fits_image(mag_file, small_file,
                                       target_resolution=None, reduction_amount=reduction)
        # print("Success!\n")
    else:
        print("Skipped! Reduced file already exists:")
        # print("\t\t", shorten_path(str(small_file), 2))
        ### WORKING HERE
        print(f"\t\tFound '{os.path.basename(small_file)}' in '\
              {shorten_path(os.path.dirname(small_file))}'")
        print("\n\t\t\t```````````````````````````````\n \n\n")

    return small_file


def reduce_fits_image(fits_path, small_file, target_resolution=None,
                      reduction_amount=None, func=np.nansum):
    """
    Open a FITS file, reduce the size of the image using astropy's block_reduce
    function, and save a new copy of the FITS file with the smaller image in the
    same directory as the original.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.
    target_resolution : int, optional
        Target resolution in pixels, if specified reduction_amount is ignored.
    reduction_amount : int, optional
        Amount to reduce the size of the image, if target_resolution is not specified.
    func : numpy function, optional
        Function to use for the reduction, defaults to np.nanmean.
    """

    print(f"\n\tReducing {fits_path}...")
    # Open the FITS file and read the data
    with fits.open(fits_path, ignore_missing_simple=True) as hdul:
        hdul.verify('silentfix')
        data = hdul[0].data
        if data is None:
            data = hdul[1].data

        current_resolution = max(data.shape)
        print("\tOriginal Shape: ", data.shape)

        # Calculate the reduction amount if target resolution is specified
        if target_resolution is not None:
            reduction_amount = int(np.ceil(current_resolution / target_resolution))

        # Raise an error if neither target_resolution nor reduction_amount is specified
        elif reduction_amount is None:
            raise ValueError("Either target_resolution or reduction_amount must be specified.")


        before_sum = np.sum(data)
        small_image = block_reduce(data, reduction_amount, func)
        after_sum = np.sum(small_image)
        if not np.isclose(before_sum, after_sum):
            print("\tREDUCTION WARNING: \n\tSum before:    ",
                  before_sum, "\n\tSum after:     ", after_sum)
        try:
            date_check =hdul[0].header["DATE"]
            useheader = hdul[0].header
        except KeyError:
            useheader = hdul[1].header

        del useheader['BLANK']
        useheader['DATAMIN'] = np.min(small_image)
        useheader['DATAMAX'] = np.max(small_image)
        useheader['BZERO'] = 0
        useheader['BSCALE'] = 1

        useheader['CDELT1'] = 360 / small_image.shape[1]  ## DEGREES
        useheader['CDELT2'] = np.deg2rad(360 / (small_image.shape[0] * np.pi)) #RADIANS

        print("\tFinal Shape:    ", small_image.shape)

        print("\tSaving  ", small_file)
        fits.writeto(small_file, small_image, useheader, overwrite=True)

        # plot_raw_magnetogram(fits_path, data, small_image)

        print("    Reduction Complete!\n")
    print("```````````````````````````\n")

    return small_file


def plot_raw_magnetogram(fits_path, data, small_image):
    """Plot the magnetogram

    Parameters
    ----------
    fits_path : str
        Magnetogram FITS file path
    data : np.ndarray
        Magnetogram data array
    small_image : np.ndarray
        Small magnetogram data array
    """

    # Save the high resolution image as a grayscale PNG
    plt.axis('off')
    high_res_output_path = fits_path.replace('.fits', '.png')
    fig = plt.gcf()
    shp = data.shape
    dmean = np.nanmean(data)
    dsig = np.nanstd(data)
    thresh = 3
    vmin = dmean - thresh*dsig
    vmax = dmean + thresh*dsig
    plt.imshow(data, cmap='gray', vmin=vmin, vmax=vmax)

    ratio = shp[1]/shp[0]
    sz0=6 #inches
    sz1=sz0*ratio #inches
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.savefig(high_res_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()

    # Save the low resolution image as a grayscale PNG
    plt.imshow(small_image, cmap='gray', vmin=vmin, vmax=vmax)
    plt.axis('off')
    low_res_output_path = fits_path.replace('.fits', '_small.png')
    fig = plt.gcf()
    shp = small_image.shape
    ratio = shp[1]/shp[0]
    sz0=6 #inches
    sz1=sz0*ratio #inches
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.savefig(low_res_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()


def load_fits_magnetogram(datdir=None, batch="fluxon", bo=2, bn=2, ret_all=False):
    """Loads a magnetogram from a FITS file.

    Parameters
    ----------
    datdir : str, optional
        Data directory path. Defaults to None.
    batch : str, optional
        Output file descriptor label. Defaults to "fluxon".
    bo : int, optional
        Output file descriptor label. Defaults to 2.
    bn : int, optional
        Output file descriptor label. Defaults to 2.
    ret_all : bool, optional
        Toggle to return both data and header. Defaults to False.

    Returns
    -------
    np.ndarray
        Output magnetogram data array
    Header
        Magnetogram header object
    """

    if not datdir:
        datdir = dat_dir
    fname = load_magnetogram_params(datdir)[2]
    fname = fname.replace("/fluxon/", f"/{batch}/").replace(f"_{bo}_", f"_{bn}_")
    fits_path = datdir + fname
    try:
        hdulist = read_fits_data(fits_path)
    except FileNotFoundError:
        hdulist = read_fits_data(fname)
    brdat = hdulist[0].data
    header= hdulist[0].header
    brdat = brdat - np.mean(brdat)
    if ret_all:
        return brdat, header
    else:
        return brdat


def write_magnetogram_params(datdir, cr, file_path, reduction):
    """Writes the magnetic_target.params file for a given CR and reduction amount.

    Parameters
    ----------
    datdir : str
        Data directory path
    cr : str
        Carrington rotation number
    file_path : str
        File path
    reduction : int
        Reduction factor
    """

    # write the parameter file
    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'w', encoding="utf-8") as fp:
        fp.write("## CR_int, Filename_str, Adapt_bool, Doplot_bool, reduction ##\n")
        fp.write(str(cr)+"\n")
        fp.write(str(file_path)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(reduction))


def load_magnetogram_params(datdir):
    """Reads the magnetic_target.params file and returns the parameters.

    Parameters
    ----------
    datdir : str
        Data directory path

    Returns
    -------
    hdr : str
        Header information
    cr : str
        Carrington rotation number
    fname : str
        Filename path
    adapt : int
        Adapt specification
    doplot : int
        Plotting toggle
    reduce : int
        Reduction factor
    """

    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'r', encoding="utf-8") as fp:
        hdr = fp.readline().rstrip()
        cr = fp.readline().rstrip()
        fname = fp.readline().rstrip()
        adapt = int(fp.readline().rstrip())
        doplot = int(fp.readline().rstrip())
        reduce = int(fp.readline().rstrip())
    return (hdr, cr, fname, adapt, doplot, reduce)


def read_fits_data(fname):
    """Reads FITS data and fixes/ignores any non-standard FITS keywords.

    Parameters
    ----------
    fname : str
        FITS file path

    Returns
    -------
    HDUList
        HDU list read from FITS file
    """
    hdulist = fits.open(fname, ignore_missing_simple=True)
    hdulist.verify('silentfix+warn')
    return hdulist


def get_fixed_coords(phi0, theta0):
    """Corrects input coordinates

    Parameters
    ----------
    phi0 : np.ndarray
        Array of phi coordinates
    theta0 : np.ndarray
        Array of theta coordinates

    Returns
    -------
    phi0: np.ndarray
        Corrected array of phi coordinates
    theta0: np.ndarray
        Corrected array of theta coordinates
    """
    ph0, th0 = phi0+np.pi, np.sin(-(theta0-(np.pi/2)))
    return ph0, th0





## Helper functions to parse the output file



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
    # print(f"WE ARE AT {this_cwd}")

    # Get the list of directories in the current directory
    dirs = list_directories(this_cwd)
    dirlist = [os.path.join(this_cwd, x) for x in dirs if "fluxon" in x]

    # Add the pipe and plotting directories to the path
    for thepath in dirlist:
        if "mhd" in thepath:
            dirlist.append(os.path.join(thepath, "flux-pipe"))
            dirlist.append(os.path.join(thepath, "flux-pipe", "plotting"))
            dirlist.append(os.path.join(thepath, "flux-pipe", "helpers"))
            break

    # Get the pipedir environment variable and add it to the path
    pipedir = os.environ.get("PIPEPATH")
    if pipedir is not None:
        dirlist.append(pipedir)

    path_add(dirlist, do_print=do_print)

    if do_print:
        print("Added fluxon directories to path.\n")

    return dirlist

def path_add(dirlist, do_print=False):    # Add the parent directory to the module search path
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