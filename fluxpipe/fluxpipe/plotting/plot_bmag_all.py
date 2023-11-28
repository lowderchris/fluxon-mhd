"""
Plotting Expansion Factor of Radial_bmag_all.dat
===============================================

This module is designed to plot the expansion factor of the given `radial_bmag_all.dat` file.
It provides a function `plot_bmag` that takes in arguments to specify the Carrington Rotation (CR),
batch name, data directory, and other parameters.

Usage:
    import plot_bmag_module
    plot_bmag_module.plot_bmag(args)

Arguments:
    args: An argparse.Namespace object containing the following attributes:
        --cr:           The Carrington Rotation for which the expansion factor is to be plotted. Default is 2163.
        --dat_dir:      The directory where the data will be stored. Default is defined in config.ini.
        --batch:        The batch name for the operation. Default is 'scalability_test'.
        --show:         Whether to show the plot or not. Default is 0.
        --nwant:        The number of fluxons wanted. Default is None.
        --file:         The file name for the data. Default is None.

Functions:
    plot_bmag(args): Plots the expansion factor based on the given arguments.

Example:
    import argparse
    args = argparse.Namespace(cr=2220, dat_dir='/path/to/data', batch='my_batch', show=1, nwant=100)
    plot_bmag_module.plot_bmag(args)

Author:
    Gilly <gilly@swri.org>

Dependencies:
    argparse, os.path, matplotlib.pyplot, numpy
"""

import argparse
import os.path
import matplotlib.pyplot as plt
import numpy as np
from fluxpipe.helpers.pipe_helper import load_fits_magnetogram, get_fixed_coords


def plot_bmag_all(args, r1=0, r2=None):
    """
    Plots the expansion factor of the given `radial_bmag_all.dat` file based on the arguments provided.

    Parameters:
        args (argparse.Namespace): An object containing the arguments for plotting.

    Returns:
        True
    """
    batch = args.batch
    filename = args.file or f'{args.dat_dir}/batches/{batch}/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    imagename = os.path.basename(filename.replace(".dat", ".png"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(filename)))
    bmagdir = os.path.join(imagedir, "imgs", "bmag")
    if not os.path.exists(bmagdir):
        os.makedirs(bmagdir)
    bmagname = os.path.join(bmagdir, imagename)

    # Load the dat file
    arr = np.loadtxt(filename).T
    nfluxon = int(arr[0,:].max())

    # Initialize Empty Arrays
    th0 = np.zeros(nfluxon)
    ph0 = np.zeros(nfluxon)
    b1  = np.zeros(nfluxon)
    th1 = np.zeros(nfluxon)
    ph1 = np.zeros(nfluxon)
    b0  = np.zeros(nfluxon)
    a0  = np.zeros(nfluxon)
    a1  = np.zeros(nfluxon)

    low_check = 1
    high_check = -4

    # Fill the arrays
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]

        rad0 = arr[4, floc[low_check]]
        rad1 = arr[4, floc[high_check]]

        th0[i] = -1 * arr[5, floc[low_check]] + np.pi/2
        th1[i] = -1 * arr[5, floc[high_check]] + np.pi/2

        ph0[i] = arr[6, floc[low_check]] + np.pi
        ph1[i] = arr[6, floc[high_check]] + np.pi

        b0[i] = arr[7, floc[low_check]]
        b1[i] = arr[7, floc[high_check]]

        a0[i] = arr[8, floc[low_check]]
        a1[i] = arr[8, floc[high_check]]

    ## Plot things!

    th0 = np.sin(th0)
    th1 = np.sin(th1)


    # Do some data manipulation
    br0_max = np.nanmax(b0) or 0.25
    br1_max = np.nanmax(b1) or 0.25
    # ar0_max = np.nanmax(a0) or 0.25
    # ar1_max = np.nanmax(a1) or 0.25

    skew = 5**2
    power = 1
    br0 = skew*(6*b0/br0_max)**power
    br1 = skew*(4*b1/br1_max)**power
    # ar0 = skew*(6*a0/ar0_max)**power
    # ar1 = skew*(4*a1/ar1_max)**power



    fig, (ax0, ax1) = plt.subplots(2)

    # Your code for the scatterplot of magnetic field values and locations


    # Your code for the line plots of radial magnetic field profiles
    RS = 696340000 #meters

    # fig, ax1 = plt.subplots()
    ax1.set_title(F"Radial Fluxon Magnetic Field Strength for CR {args.cr}")

    alabel = r"Magnetic Field Strength [Gauss]"
    blabel = r"Expansion [Ratio]"


    ax1.set_xlabel(r'Heliocentric Radius [R$_\odot$]')
    ax1.set_ylabel(alabel, color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    # ax1.set_yscale('log')
    ax1.set_xscale('log')

    # Create a second y-axis for the expansion data
    ax2 = ax1.twinx()
    # ax2.set_yscale('log')
    ax2.set_ylabel(blabel, color='b')
    ax2.tick_params(axis='y', labelcolor='b')

    from matplotlib.colors import LinearSegmentedColormap

    # Creating custom colormaps
    # colors_red = [(1, 0, 0, i/nfluxon) for i in range(nfluxon)]  # Dark to light red
    # colors_blue = [(0, 0, 1, i/nfluxon) for i in range(nfluxon)]  # Dark to light blue
    # custom_cmap_red = LinearSegmentedColormap.from_list("custom_red", colors_red, N=nfluxon)
    # custom_cmap_blue = LinearSegmentedColormap.from_list("custom_blue", colors_blue, N=nfluxon)

    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]
        rr = arr[4,floc]/RS
        field = arr[7,floc]
        expansion = arr[8,floc]

        # Use the 'viridis' colormap for the field plot
        ax1.plot(rr, field, c=plt.cm.Reds(i / nfluxon), label=blabel if i == 0 else "", alpha=0.8)

        # Use the 'plasma' colormap for the expansion plot
        ax2.plot(rr, expansion, c=plt.cm.Blues(i / nfluxon), label=alabel if i == 0 else "", alpha=0.8)

    # You might need to adjust the scale for ax2 based on the range of your expansion data
    # ax2.set_yscale('log')

    r1 = rad0/RS
    r2 = rad1/RS
    ax1.axvline(r1, ls=":", c='k', zorder=1000)
    ax1.axvline(r2, ls=":", c='k', zorder=1000)




    magnet = load_fits_magnetogram(configs=configs)
    ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto')
    sc00 = ax0.scatter(ph0, th0, c=b0, s = b0, cmap="winter",
                    alpha=0.75, label=r"B(1.0Rs)")
    sc01 = ax0.scatter(ph1, th1, c=b1, s = 10*b1, cmap="autumn",
                    alpha=0.75, label=r"B(21.5Rs)", marker='s')
    cbar01 = fig.colorbar(sc01, ax=ax0)
    cbar00 = fig.colorbar(sc00, ax=ax0)
    cbar00.set_label(r"B(1.0Rs)  [Gauss]")
    cbar01.set_label(r"B(21.5Rs) [Gauss]")

    fig.suptitle(F"Fluxon Magnetic Field Strength for CR {args.cr}")
    ax0.set_title(F"Magnetic Field Strength by Location")



    # Adjust layout
    fig.set_size_inches((8,8))
    plt.tight_layout()
    plt.savefig(bmagname)
    if args.show or True:
        plt.show()
    plt.close(fig)
    print("Done!\n")
    return True

# Main Code
if __name__ == "__main__":
    # Create the argument parser
    print("\tPlotting Bmag_All...", end="")
    from fluxpipe.helpers.pipe_helper import configurations
    configs = configurations()

    parser = argparse.ArgumentParser(description=
                            'This script plots the expansion factor of the given radial_bmag_all.dat')
    parser.add_argument('--cr',     type=int, default=configs['rotations'][0],    help='Carrington Rotation')
    parser.add_argument('--dat_dir',type=str, default=configs["data_dir"],        help='data directory')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"],      help='select the batch name')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='magnetogram file')
    parser.add_argument('--show',   type=int, default=0)
    parser.add_argument('--file',   type=str, default=None)
    parser.add_argument('--adapt',  type=int, default=configs["adapt"],           help='Use ADAPT magnetograms')

    args = parser.parse_args()

    plot_bmag_all(args)
