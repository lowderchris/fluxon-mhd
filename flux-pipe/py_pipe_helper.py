"""
py_pipe_helper: Library of Helper Scripts for the Flux-Pipe Algorithm
=====================================================================

This library contains a collection of utility functions designed to assist with
the flux-pipe algorithm. It includes methods for managing directories, handling
FITS files, and manipulating magnetogram data.

Usage:
    import flux_pipe_helper as fph

Arguments:
    No command-line arguments are required when using this library as a module.

Functions:
    add_dir_to_path, add_top_level_dirs_to_path, add_paths, make_mag_dir, get_magnetogram_file,
    reduce_mag_file, reduce_fits_image, plot_raw_magnetogram, load_fits_magnetogram,
    write_magnetogram_params, load_magnetogram_params, find_file_with_string, shorten_path,
    read_fits_data, get_fixed_coords

Example:
    import flux_pipe_helper as fph
    fph.add_dir_to_path("/path/to/directory")

Author:
    Gilly <gilly@swri.org> (and others!)

Dependencies:
    os, os.path, sys, pathlib.PosixPath, Path, plotting.py_plot_helper, numpy as np,
    matplotlib.pyplot as plt, astropy.io.fits, astropy.nddata.block_reduce, sunpy,
    sunpy.coordinates, sunpy.io, sunpy.net.Fido, attrs as a
"""


# Import libraries
import os
import os.path
import sys
from pathlib import PosixPath, Path
import plotting.py_plot_helper

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.nddata import block_reduce

import sunpy
import sunpy.coordinates
import sunpy.io
from sunpy.net import Fido, attrs as a

default_email = "chris.gilly@colorado.edu"
default_root_dir = "/Users/cgilbert/vscode/fluxons/fluxon-mhd/"
dat_dir = "/Users/cgilbert/vscode/fluxons/fluxon-data/"

def add_dir_to_path(root_dir=None):
    """Adds a directory and all subdirectories to the PATH environment variable.

    Parameters
    ----------
    root_dir : str, optional
        Root directory path
    """

    if root_dir is None:
        root_dir = default_root_dir

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
