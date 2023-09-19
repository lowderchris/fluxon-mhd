"""
Plotting Expansion Factor of Radial_fr.dat
==========================================

This module is designed to plot the expansion factor of the given `radial_fr.dat` file.
It provides a function `plot_fr` that takes in arguments to specify the Carrington Rotation (CR),
batch name, data directory, and other parameters.

Usage:
    import plot_fr_module
    plot_fr_module.plot_fr(args)

Arguments:
    args: An argparse.Namespace object containing the following attributes:
        --cr:           The Carrington Rotation for which the expansion factor is to be plotted. Default is 2163.
        --dat_dir:      The directory where the data will be stored. Default is defined in config.ini.
        --batch:        The batch name for the operation. Default is 'scalability_test'.
        --show:         Whether to show the plot or not. Default is 0.
        --nwant:        The number of fluxons wanted. Default is None.
        --file:         The file name for the data. Default is None.

Functions:
    plot_fr(args): Plots the expansion factor based on the given arguments.

Example:
    import argparse
    args = argparse.Namespace(cr=2220, dat_dir='/path/to/data', batch='my_batch', show=1, nwant=100)
    plot_fr_module.plot_fr(args)

Author:
    Gilly <gilly@swri.org>

Dependencies:
    argparse, os.path, matplotlib.pyplot, numpy
"""

import argparse
import os.path
import matplotlib.pyplot as plt
import numpy as np


def plot_fr(args):
    """
    Plots the expansion factor of the given `radial_fr.dat` file based on the arguments provided.

    Parameters:
        args (argparse.Namespace): An object containing the arguments for plotting.

    Returns:
        True
    """
    batch = args.batch
    filename = args.file or f'{args.dat_dir}/batches/{batch}/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_fr.dat'
    imagename = os.path.basename(filename.replace(".dat", ".png"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(filename)))
    frdir = os.path.join(imagedir, "imgs", "fr")
    if not os.path.exists(frdir):
        os.makedirs(frdir)
    frname = os.path.join(frdir, imagename)

    # Load the dat file
    arr = np.loadtxt(filename).T
    nfluxon = int(arr[0,:].max())

    # Initialize Empty Arrays
    th0 = np.zeros(nfluxon)
    ph0 = np.zeros(nfluxon)
    a0  = np.zeros(nfluxon)
    th1 = np.zeros(nfluxon)
    ph1 = np.zeros(nfluxon)
    a1  = np.zeros(nfluxon)

    # Fill the arrays
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]

        th0[i] = -1 * arr[5, floc[0]] + np.pi/2
        th1[i] = -1 * arr[5, floc[-1]] + np.pi/2

        ph0[i] = arr[6, floc[0]] + np.pi
        ph1[i] = arr[6, floc[-1]] + np.pi

        a0[i] = arr[7, floc[0]]
        a1[i] = arr[7, floc[-1]]

    ## Plot things!

    fig, (ax1, ax2) = plt.subplots(2)
    # TODO make sure the X axis is correct
    shift = -1.0
    ax1.set_title(F"Fluxon Expansion Factors for CR {args.cr}")
    RS = 696340000 #meters
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]
        ax1.plot((shift + arr[4,floc]/RS), arr[7,floc])
    r1 = 1.
    r2 = 21.5
    ax1.axvline(r1, ls=":", c='lightgrey')
    ax1.axvline(r2, ls=":", c='lightgrey')
    ax1.set_xlabel(r'Heliocentric Radius [R$_\odot$]')
    ax1.set_ylabel('Cross Sectional Area [m$^2$]')
    ax1.set_yscale('log')
    ax1.set_ylim((10.0**10, 10.0**20))

    ax2.set_title(r"Fluxon Area at Low (1R$_\odot$) and High (21.5R$_\odot$) Radii")

    a0_max = np.nanmax(np.abs(a0)) or 0.25
    a1_max = np.nanmax(np.abs(a1)) or 0.25
    skew = 2

    aa0 = (skew+5*np.abs(a0)/a0_max)**3
    aa1 = (skew+5*np.abs(a1)/a1_max)**3

    ax2.scatter(ph0, np.sin(th0), s = aa0 , label=r'A(1.0R$_\odot$),', alpha=0.75, marker='s')
    ax2.scatter(ph1, np.sin(th1), s = aa1 , label=r'A(21.5R$_\odot$)', alpha=0.75, marker='s')
    ax2.set_xlabel('Longitude (Radians)')
    ax2.set_ylabel('Sine latitude')
    ax2.set_ylim((-1.1,1.1))
    ax2.axhline(-1, c='lightgrey', zorder=-10)
    ax2.axhline( 1, c='lightgrey', zorder=-10)
    ax2.axvline(0, c='lightgrey', zorder=-10)
    ax2.axvline(2*np.pi, c='lightgrey', zorder=-10)
    ax2.legend(loc="upper right")

    fig.set_size_inches((6,8))
    plt.tight_layout()
    plt.savefig(frname)
    if args.show:
        plt.show()
    plt.close(fig)
    print("Done!")
    return True




########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    # Create the argument parser
    print("\tPlotting Fr...", end="")
    from pipe_helper import configurations
    configs = configurations()

    parser = argparse.ArgumentParser(description=
                            'This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr',     type=int, default=configs['rotations'][0],    help='Carrington Rotation')
    parser.add_argument('--dat_dir',type=str, default=configs["data_dir"],        help='data directory')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"],      help='select the batch name')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='magnetogram file')
    parser.add_argument('--show',   type=int, default=0)
    parser.add_argument('--file',   type=str, default=None)

    args = parser.parse_args()

    plot_fr(args)