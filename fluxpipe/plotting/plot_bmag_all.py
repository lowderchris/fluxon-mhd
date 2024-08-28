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
from matplotlib import gridspec
import numpy as np
import shutil
from fluxpipe.helpers.pipe_helper import load_fits_magnetogram, get_fixed_coords
import sunpy.coordinates
import cv2


def plot_bmag_all(args, r1=1, r2=-2, do_r2=False, maxlist=None):
    """
    Plots the expansion factor of the given `radial_bmag_all.dat` file based on the arguments provided.

    Parameters:
        args (argparse.Namespace): An object containing the arguments for plotting.

    Returns:
        True
    """
    batch = args.batch
    filename = args.file or f'{args.dat_dir}/batches/{batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    imagename = os.path.basename(filename.replace(".dat", f"_{r1:02d}.png"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(filename))))
    bmagdir = os.path.join(imagedir, "imgs", "bmag", "all", "sets", f"cr{args.cr}_all_f{args.nwant}")
    if not os.path.exists(bmagdir):
        os.makedirs(bmagdir)
    bmagname = os.path.join(bmagdir, imagename)


    ### GET TOP AND BOTTOM DATA
    # Load the dat file
    arr = np.loadtxt(filename, skiprows=1).T
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

    low_check = r1
    high_check = r2

    # Fill the arrays
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]
        lower = floc[low_check]
        upper = floc[high_check]

        rad0 = arr[4, lower]
        rad1 = arr[4, upper]

        th0[i] = -1 * arr[5, lower] + np.pi/2
        th1[i] = -1 * arr[5, upper] + np.pi/2

        ph0[i] = np.mod(arr[6, lower], 2*np.pi)
        ph1[i] = np.mod(arr[6, upper], 2*np.pi)

        a0[i] = arr[7, lower]
        a1[i] = arr[7, upper]

        b0[i] = arr[8, lower]
        b1[i] = arr[8, upper]

        # b0x[i] = arr[9, lower]
        # b0y[i] = arr[10, lower]
        # b0z[i] = arr[11, lower]

        # b1x[i] = arr[9, upper]
        # b1y[i] = arr[10, upper]
        # b1z[i] = arr[11, upper]

    ## Plot things!

    th0 = np.sin(th0)
    th1 = np.sin(th1)

    if maxlist is None:
        maxlist = th0

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







    # fig, (ax0, ax1) = plt.subplots(2)

    RS = 696340000 #meters


    # Using GridSpec for layout
    fig = plt.figure(figsize=(10, 8))  # Adjust the figure size as needed
    gs = gridspec.GridSpec(3, 1, height_ratios=[4, 4, 1])  # Adjust the ratios as per your requirement

    # Magnetogram plot
    ax0 = plt.subplot(gs[0])

    # Magnetic Field Strength and Expansion Factor plot
    ax1 = plt.subplot(gs[1])
    ax2 = ax1.twinx()

    # Sunspot Number plot
    ax3 = plt.subplot(gs[2])




    ### FIRST PLOT ###
    # Plot magnetogram
    magnet = load_fits_magnetogram(configs=configs)
    ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=0.5)


    #Plot Magnetic Field Strength
    try:
        bb0 = (b0 - np.nanmin(b0)) / (np.nanmax(b0) - np.nanmin(b0))
        if np.isnan(bb0).any():
            bb0[np.isnan(bb0)] = b0[np.isnan(bb0)] / np.nanmax(b0)
    except (ValueError, RuntimeWarning):
        bb0 = b0 / np.nanmax(b0)

    # bb0 = (bb0 + 0.2) / 1.2

    if r1 and False:
        # sc00 = ax0.scatter(ph0, th0, c=b0, s=50, cmap="binary_r", alpha=0.8,
        sc00 = ax0.scatter(ph0, th0, c="k", s=50, alpha=1-bb0,
                           label=r"B(1.0Rs)", zorder=100) #, vmin=0, vmax=25)
        cbar00 = fig.colorbar(sc00, ax=ax0)
        cbar00.set_label(f"B({r1:.2f} Rs)  [Gauss]")
        pass

    # Plot Latitude of points
    do_lat = True
    if do_lat:
        sc01 = ax0.scatter(ph0, th0, c=maxlist, s = 100, cmap="brg",
                    alpha=bb0, label=r"First Point Latitude", zorder = 99, marker='o', vmin=-1, vmax=1)
        cbar01 = fig.colorbar(sc01, ax=ax0)
        cbar01.set_label(f"First Point Sine Latitude")

    if do_r2:
        sc01 = ax0.scatter(ph1, th1, c=b1, s = 10*b1, cmap="autumn",
                    alpha=0.75, label=r"B(21.5Rs)", marker='s', vmin=0, vmax=25)
        cbar01 = fig.colorbar(sc01, ax=ax0)
        # cbar01.set_label(f"B({r2:.2f} Rs) [Gauss]")
        cbar01.set_label("First Point Sine Latitude")







    ### SECOND PLOT ###
    alabel = r"Magnetic Field Strength [Gauss]"
    blabel = r"Expansion Factor"
    tlabel = r"Latitude [Radians]"


    ax1.set_xlabel(r'Heliocentric Radius [R$_\odot$]')
    ax1.set_ylabel(alabel, color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    ax1.set_yscale('log')
    ax1.set_xscale('log')

    # Create a second y-axis for the expansion data
    ax2.set_yscale('log')
    ax2.set_ylabel(blabel, color='b')
    ax2.tick_params(axis='y', labelcolor='b')


    # Plot the Curves
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]
        rr = arr[4, floc]

        phi = arr[6, floc]
        ttheta = -1 * arr[5, floc] + np.pi/2
        pphi = np.mod(arr[6, floc], 2*np.pi)


        expansion = arr[7,floc]
        field = arr[8,floc]
        field_2 = arr[9,floc]

        ax1.plot(rr-1, field, c=plt.cm.Reds(i / nfluxon), label=blabel if i == 0 else "", alpha=0.8)
        ax2.plot(rr-1, expansion, c=plt.cm.Blues(i / nfluxon), label=alabel if i == 0 else "", alpha=0.8)
        # ax1.plot(rr, pphi, c=plt.cm.Purples(i / nfluxon), label=tlabel if i == 0 else "", alpha=0.8)

    r1 = rad0
    r2 = rad1
    ax1.axvline(r1-1, ls=":", c='k', zorder=1000000)
    if do_r2:
        ax1.axvline(r2-1, ls=":", c='k', zorder=1000)





    ### THIRD PLOT ###
    # Plot the Sunspot Number
    carrington = np.loadtxt("/Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/SN_m_tot_V2.0.txt").T
    ## https://sidc.be/SILSO/datafiles#total ##
    import sunpy.coordinates
    date = carrington[2]
    sunspots = carrington[3]
    this_date = sunpy.coordinates.sun.carrington_rotation_time(args.cr)
    # CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))

    # fig, ax = plt.subplots()
    ax3.plot(date, sunspots, label="Sunspots")
    ax3.axvline(this_date.decimalyear, ls=":", c='k', zorder=1000000)
    ax3.set_xlabel("Year")
    ax3.set_ylabel("Sunspots")
    ax3.set_title("Solar Cycle")
    ax3.set_xlim(2005, 2025)
    ax3.set_ylim(0, 200)
    # fig.close()
    # plt.show()





    fig.suptitle(F"Fluxon Magnetic Field Properties for CR {args.cr}, with {nfluxon} Open Fields")
    ax0.set_title(F"Field Strength and Origin at r = {r1:.2f} Rs")
    ax1.set_title(F"Magnetic Field Strength and Expansion Factor vs Height")

    ax0.set_xlabel('Longitude (Radians)')
    ax0.set_ylabel('Sine latitude')

    # Adjust layout
    fig.set_size_inches((8,8))
    plt.tight_layout()
    plt.savefig(bmagname)
    if args.show or True:
        plt.show()
    plt.close(fig)
    print(".", end="", flush=True)
    return bmagname, maxlist



import cv2
import os

def create_video_from_images(image_folder, video_name, frame_rate=3.0, codec='XVID', ttype="png"):
    """
    Create a video from a sequence of images in a folder.

    :param image_folder: Path to the directory containing images.
    :param video_name: Path where the video will be saved.
    :param frame_rate: Frame rate of the output video.
    :param codec: Codec to be used for creating the video.
    """
    images = [img for img in os.listdir(image_folder) if img.endswith(ttype)]
    if not images:
        raise ValueError("No images found in the folder.")

    images.sort()  # Sort the images

    # Determine the width and height from the first image
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    if frame is None:
        raise ValueError("Unable to read the first image.")

    height, width, layers = frame.shape

    video_folder = os.path.join(os.path.dirname(image_folder), "expansion_videos")
    if not os.path.exists(video_folder):
        os.makedirs(video_folder)

    # Initialize video writer
    fourcc = cv2.VideoWriter_fourcc(*codec)
    video = cv2.VideoWriter(os.path.join(video_folder, video_name), fourcc, frame_rate, (width, height))

    # Add images to video
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    video.release()


def run_plots(args, times=0):
    maxlist = None
    for rr in np.arange(times+1):
        tried = True
        try:
            bmagfull, maxlist = plot_bmag_all(args, rr, maxlist=maxlist)
            bmagdir = os.path.dirname(bmagfull)
            bmagfile= os.path.basename(bmagfull)
            if rr == 0:
                new_dir = os.path.join(bmagdir, "../..", "bottom")
                if not os.path.exists(new_dir):
                    os.makedirs(new_dir)
                shutil.copyfile(bmagfull, os.path.join(new_dir, bmagfile))
        except IndexError:
            new_dir = os.path.join(bmagdir, "../..", "top")
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            shutil.copyfile(bmagfull, os.path.join(new_dir, bmagfile))
            break
    assert tried, "Something went wrong with the plot_bmag_all function."
    print("Done!\n")
    return bmagdir


# Main Code
if __name__ == "__main__":
    # Create the argument parser
    print("\n\tPlotting Bmag_All...", end="")
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


    bmagdir = run_plots(args)
    filename = "expansion_cr{}_f{}.avi".format(args.cr, args.nwant)
    # create_video_from_images(bmagdir, filename)

