"""
Plotting Magnetic Field Strength and Fluxon Area
================================================

This script plots the magnetic field strength and fluxon area of the fluxons
at the lower and upper boundaries. It provides options to specify the Carrington Rotation (CR),
batch name, data directory, and other parameters.

Usage:
    python plot_bmag.py [--cr CARRINGTON_ROTATION]
                        [--show SHOW_PLOT] [--batch BATCH_NAME]
                        [--file FILE_NAME] [--nwant NUMBER_OF_FLUXONS]

Parameters:
    --cr:       Carrington Rotation for which the data is to be plotted. Default is 0.
    --show:     Boolean flag to indicate whether to show the plot or not. Default is 0.
    --batch:    Batch name for the operation. Default is 'default_batch'.
    --file:     File name for the data. Default is a constructed path based on other parameters.
    --nwant:    Number of fluxons to plot. Default is None.

Functions:
    None (script-based)

Example:
    python plot_bmag.py --cr 2183 --show 1 --batch 'my_batch' --nwant 100

Dependencies:
    os.path, argparse, matplotlib.pyplot, numpy, pipe_helper

Author:
    Gilly <gilly@swri.org> (and others!)

"""


import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pipe_helper import (configurations, get_fixed_coords, load_fits_magnetogram,
                         load_magnetogram_params, shorten_path)


def plot_bmag(configs):

    filename = configs.get('file', f"{configs['data_dir']}/batches/{configs['batch_name']}/cr{configs['cr']}/wind/cr{configs['cr']}_f{configs['nwant']}_radial_bmag.dat")

    if configs['adapt'] and not "adapt" in filename:
        filename = filename.replace(configs['batch_name'], f"{configs['batch_name']}_adapt")


    # Load the dat file
    arr = np.loadtxt(filename).T
    fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1 = arr

    # Convert coords to correct coords
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)

    # Do some data manipulation
    br0_max = np.nanmax(br0) or 0.25
    br1_max = np.nanmax(br1) or 0.25
    ar0_max = np.nanmax(ar0) or 0.25
    ar1_max = np.nanmax(ar1) or 0.25

    skew = 5**2
    power = 1
    b0 = skew*(6*br0/br0_max)**power
    b1 = skew*(4*br1/br1_max)**power
    a0 = skew*(6*ar0/ar0_max)**power
    a1 = skew*(4*ar1/ar1_max)**power

    # Plot the Data
    fig, (ax0, ax1) = plt.subplots(2)

    magnet = load_fits_magnetogram(cr=configs['cr'], adapt=configs['adapt'],    )
    ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto')
    ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto')

    sc00 = ax0.scatter(ph0, th0, c=br0, s = b0, cmap="winter",
                    alpha=0.75, label=r'B(1.0R$_\{odot}$)')
    sc01 = ax0.scatter(ph1, th1, c=br1, s = b1, cmap="autumn",
                    alpha=0.75, label=r'B(21.5R$_\{odot}$)', marker='s')
    sc10 = ax1.scatter(ph0, th0, c=ar0, s = a0, cmap="winter",
                    alpha=0.75, label=r'A(1.0R$_\{odot}$)')
    sc11 = ax1.scatter(ph1, th1, c=ar1, s = a1, cmap="autumn",
                    alpha=0.75, label=r'A(21.5R$_\{odot}$)', marker='s')

    cbar01 = fig.colorbar(sc01, ax=ax0)
    cbar00 = fig.colorbar(sc00, ax=ax0)
    cbar11 = fig.colorbar(sc11, ax=ax1)
    cbar10 = fig.colorbar(sc10, ax=ax1)

    cbar00.set_label(r"B(1.0R$_\{odot}$)  [Gauss] ")
    cbar01.set_label(r"B(21.5R$_\{odot}$) [Gauss]")
    cbar10.set_label(r"A(1.0R$_\{odot}$)  [m$^2$]")
    cbar11.set_label(r"A(21.5R$_\{odot}$) [m$^2$]")

    ax0.set_title(F"Fluxon Magnetic Field Strength for CR {args.cr}")
    ax1.set_title(F"Fluxon Area for CR {args.cr}")

    for ax in (ax0, ax1):
        ax.set_xlabel('Longitude (Radians)')
        ax.set_ylabel('Sine latitude')
        ax.set_ylim((-1.5,1.1))
        ax.axhline(-1, c='lightgrey', zorder=-10)
        ax.axhline( 1, c='lightgrey', zorder=-10)
        ax.axvline(0, c='lightgrey', zorder=-10)
        ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        ax.legend(loc="lower right")

    fig.set_size_inches((8,8))
    plt.tight_layout()

    imagename = os.path.basename(filename.replace(".dat", ".png"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(filename)))
    bdir = os.path.join(imagedir, "imgs", "bmag")
    if not os.path.exists(bdir):
        os.makedirs(bdir)
    pngname = os.path.join(bdir, imagename)

    plt.savefig(pngname)
    if configs['verbose']:
        plt.show()
    plt.close(fig)
    print("Done!\n")
    print("\t\tSaved to", shorten_path(pngname), "\n")




########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    # Create the argument parser
    print("\n\tPlotting Bmag...", end="")
    parser = argparse.ArgumentParser(description=
                                    'This script plots the radial magnetic field')
    parser.add_argument('--cr', type=int, default=2219, help='Carrington Rotation')
    parser.add_argument('--file', type=str, default=None, help='Data File Name')
    parser.add_argument('--nwant', type=int, default=100, help='Number of Fluxons')
    parser.add_argument('--adapt', type=int, default=1, help='using adapt maps')
    args = parser.parse_args()
    configs = configurations(debug=False, args=args)
    plot_bmag(configs)
