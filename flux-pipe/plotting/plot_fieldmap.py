"""
Primary Script for Plotting the Magnetogram with Footpoints
==========================================================

This script is designed to plot the magnetogram along with footpoints. It accepts
configuration parameters through a dictionary and command-line arguments through an argparse.Namespace object.

Usage:
    python plot_fieldmap.py [--cr CARRINGTON_ROTATION] [--nwant NUMBER_WANTED]
                             [--open OPEN_FLUXONS] [--closed CLOSED_FLUXONS]

Arguments:
    --cr, --nwant, --open, --closed, --dat_dir, --batch:
        These can now be passed through the 'config' dictionary and 'args' argparse.Namespace object.

Functions:
    magnet_plot:    The primary function for plotting the magnetogram with footpoints.

Example:
    python plot_fieldmap.py --cr 2220 --nwant 100 --open 50 --closed 50 --dat_dir '/path/to/data' --batch 'my_batch'

Author:
    Gilly <gilly@swri.org> (and others!)

Dependencies:
    os, os.path, argparse, matplotlib.pyplot, numpy, py_plot_helper, pfss_funcs, py_pipe_helper
"""

import os
import os.path as path
import argparse
import matplotlib.pyplot as plt
import numpy as np

from py_plot_helper import get_ax
from pfss_funcs import pixel_to_latlon
from py_pipe_helper import (load_fits_magnetogram, load_magnetogram_params, shorten_path)

def magnet_plot(config, args, force=False, nact=0, do_print_top=False, ax=None, verb=True, ext="pdf",
                plot_all=True, plot_open=True, do_print=False, vmin=-500, vmax=500):
    """
    The primary function for plotting the magnetogram with footpoints.

    Parameters
    ----------
    config : dict
        Dictionary containing configuration parameters.
        - CR : int, the Carrington rotation number
        - data_dir : str, the data directory
        - batch_name : str, the name of the batch
        - mag_reduct : int, optional, the factor by which to reduce the magnetogram, default is 0
        - nwant : int, optional, the desired number of fluxons, default is None

    args : argparse.Namespace
        Argument parser namespace containing command-line arguments.
        - open : int, optional, number of open fluxons, default is None
        - closed : int, optional, number of closed fluxons, default is None

    force : bool, optional
        Force the plot to occur, default is False
    nact : int, optional
        The actual number of fluxons, default is 0
    do_print_top : bool, optional
        Print the top, default is False
    ax : pyplot axis, optional
        The axis to plot upon, default is None
    verb : bool, optional
        Whether to be verbose in output, default is True
    ext : str, optional
        The type of plot to save, default is "pdf"
    plot_all : bool, optional
        Whether to plot all the fluxons, default is True
    plot_open : bool, optional
        Whether to plot the open fields, default is True
    do_print : bool, optional
        Even more verbosity, default is False
    vmin : int, optional
        Minimum value of the colortable, default is -500
    vmax : int, optional
        Maximum value of the colortable, default is 500

    Returns
    -------
    _n_open, _n_closed, _n_flux, _fnum, _n_outliers
        Integers describing the number of open fluxons, closed fluxons, total fluxons,
        actual fluxons, and outliers, respectively
    """
    # Unpack configuration and args parameters
    datdir = config["data_dir"]
    _batch = config["batch_name"]
    reduce_amt = config["mag_reduct"]
    get_cr = args.cr
    nwant = args.nwant
    open_f = args.open
    closed_f = args.closed

    fig, ax0 = get_ax(ax)
    if do_print:
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("(py) Plotting Magnetic Field Maps")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


    if do_print_top:
        print("\n\tMaking Magnetogram with Footpoints...")
    # Define the directory paths for the files
    floc_path = f"{datdir}/batches/{_batch}/cr{get_cr}/floc/"
    top_dir   = f"{datdir}/batches/{_batch}/imgs/footpoints/"
    if not path.exists(top_dir):
        os.makedirs(top_dir)

    # Define the file names with their complete paths
    open_file   = open_f     or   f"{floc_path}floc_open_cr{get_cr}_r{reduce_amt}_f{nwant}.dat"
    closed_file = closed_f or     f"{floc_path}floc_closed_cr{get_cr}_r{reduce_amt}_f{nwant}.dat"
    all_file    = closed_file.replace("closed_", "")

    # Load the data
    if do_print_top:
        print(f"\t\tOpening {shorten_path(all_file)}...")
    fluxon_location = np.genfromtxt(all_file)

    if do_print_top:
        print(f"\t\tOpening {shorten_path(open_file)}...")
    oflnum, oflx, olat, olon, orad = np.loadtxt(open_file, unpack=True)

    if do_print_top:
        print(f"\t\tOpening {shorten_path(closed_file)}...\n")
    cflnum, _, _, _, crad = np.loadtxt(closed_file, unpack=True)

    magnet, header = load_fits_magnetogram(batch=_batch, ret_all=True)
    f_lat, f_lon, f_sgn, _fnum = pixel_to_latlon(magnet, header, fluxon_location)

    ## Keep only the values where the radius is 1.0
    rtol = 0.01
    get_r = 1.0

    #Open fields
    oflnum_low = oflnum[np.isclose(orad, get_r, rtol)]
    oflx_low =     oflx[np.isclose(orad, get_r, rtol)]
    olat_low =     olat[np.isclose(orad, get_r, rtol)]
    olon_low =     olon[np.isclose(orad, get_r, rtol)]

    # Closed fields
    cflnum_low = cflnum[np.isclose(crad, get_r, rtol)]


    # Convert to radians
    ph_olow, th_olow = np.sin(np.deg2rad(olat_low)), np.deg2rad(olon_low)

    # Report the number of open and closed fluxons
    _n_open = int(np.max(oflnum_low))
    _n_closed = int(np.max(cflnum_low))
    _n_flux = _n_open + _n_closed
    _n_outliers = np.abs(_fnum-_n_flux)

    # Define the file name for the plot
    pic_name = f'cr{get_cr}_f{nwant}_ou{_n_open}_footpoints_topology.{ext}'
    fluxon_map_output_path =   path.join(floc_path, pic_name)
    fluxon_map_output_path_top = path.join(top_dir, pic_name)

    # Check if the plot already exists
    do_plot = False
    pic_paths = [fluxon_map_output_path, fluxon_map_output_path_top]
    # pic_paths = [fluxon_map_output_path_top]
    for testpath in pic_paths:
        if not path.exists(testpath):
            do_plot = True
            break

    if do_print:
        print("\tPlotting...", end="")
    if do_plot or force or (ax is not None):
        # Plot the magnetogram

        ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
                extent=(0,2*np.pi,-1,1), aspect='auto', vmin=vmin, vmax=vmax)

        # Plot all the fluxons
        if plot_all:
            colors_all = ['orange' if s > 0 else 'teal' for s in f_sgn]
            ax0.scatter(f_lon, f_lat, s=3**2, c=colors_all, alpha=0.6)

        # Plot the open fluxons
        if plot_open:
            colors_open = ['red' if s > 0 else 'blue' for s in oflx_low]
            ax0.scatter(th_olow, ph_olow, s=5**2, c=colors_open,
                               alpha=1.0, label='Open Fields', edgecolors='k')

        if ax is None:
            shp = magnet.shape #pixels
            plt.axis('off')
            sz0=6 #inches
            ratio = shp[1]/shp[0]
            sz1=sz0*ratio #inches
            DPI = shp[1] / sz1 #pixels/inch
            fig.set_size_inches((sz1, sz0))
            plt.tight_layout()
            plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)
            # plt.show()
            plt.close(fig)
    else:
        if do_print:
            print("\tSkipped! Files already exist:")
            print(f"\t\t{shorten_path(fluxon_map_output_path)}")
            print(f"\t\t{shorten_path(fluxon_map_output_path_top)}")
    if do_print:
        print(f"\n\t    n_open: {_n_open}, n_closed: {_n_closed}, \
                n_total: {_n_flux}, n_all: {_fnum}, n_outliers: {_n_outliers}")

    if do_print_top:
        print("\t\t    Success!")
        print("\t\t\t```````````````````````````````\n\n")
    return _n_open, _n_closed, _n_flux, _fnum, _n_outliers





########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
                                    'This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--nwant', type=int, default=None, help='Number of fluxons to create')
    parser.add_argument('--file', type=str, default=None, help='path to the dat file')
    parser.add_argument('--open', type=str, default=None)
    parser.add_argument('--closed', type=str, default=None)
    args = parser.parse_args()

    from config_reader import load_configs

    configs = load_configs()

    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(
    config=configs,
    args=args,
    do_print=True,
    do_print_top=True
)


    # n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(configs["CR"], configs["data_dir"], configs["batch_name"], args.open,
    #         args.closed, do_print=True, reduce_amt=configs["mag_reduct"], nwant=configs["nwant"], do_print_top=True)
