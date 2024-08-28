"""
!!NEEDS UPDATE!!: Plotting Expansion Factor of Radial_bmag_all.dat
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

def_file = "/Users/cgilbert/vscode/fluxons/fluxon-data/zephyr_2007_2013.sav"

def load_zephyr(file = def_file):
    from scipy.io import readsav
    data = readsav(file)
    # nz ()
    # nmods ()
    # model_year (319,)
    # tag1 (319,)
    # tag2 (319,)
    # rx (1300,)
    # rho (319, 1300)
    # uu (319, 1300)
    # valf (319, 1300)
    # t (319, 1300)
    # br (319, 1300)
    # [print(k, data[k].shape) for k in data.keys()]

    return data


def sunspotplot(ax3, crlist=None):
    ### THIRD PLOT ###
    # Plot the Sunspot Number
    carrington = np.loadtxt("/Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/SN_m_tot_V2.0.txt").T
    ## https://sidc.be/SILSO/datafiles#total ##
    import sunpy.coordinates
    date = carrington[2]
    sunspots = carrington[3]
    this_date = sunpy.coordinates.sun.carrington_rotation_time(args.cr)
    ax3.axvline(this_date.decimalyear, ls=":", c='k', zorder=1000000)

    if crlist:
        for cr in crlist:
            this_date = sunpy.coordinates.sun.carrington_rotation_time(cr)
            ax3.axvline(this_date.decimalyear, ls=":", c='k', zorder=1000000)

    # CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))

    ax3.plot(date, sunspots, label="Sunspots")
    ax3.set_xlabel("Year")
    ax3.set_ylabel("Sunspots")
    ax3.set_title("Solar Cycle Phase")
    ax3.set_xlim(2005, 2025)
    ax3.set_ylim(0, 200)

def multiplot_bmag_all(args, r1=1, r2=-1, do_r2=False, maxlist=None, do_stiff=False):
    """
    Plots the expansion factor of the given `radial_bmag_all.dat` file based on the arguments provided.

    Parameters:
        args (argparse.Namespace): An object containing the arguments for plotting.

    Returns:
        True
    """

    # Parse the arguments and directories
    batch = args.batch
    relaxed_filename = args.file or f'{args.dat_dir}/batches/{batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    stiff_filename = f'{args.dat_dir}/batches/{batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_stiff_radial_bmag_all.dat'

    if do_stiff:
        filename = stiff_filename
    else:
        filename = relaxed_filename

    imagename = os.path.basename(filename.replace(".dat", f"_{r1:02d}.png"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(filename))))
    bmagdir = os.path.join(imagedir, "imgs", "bmag", "all", "sets", f"cr{args.cr}_all_f{args.nwant}")
    if not os.path.exists(bmagdir):
        os.makedirs(bmagdir)
    bmagname = os.path.join(bmagdir, imagename)


    # Using GridSpec for layout
    fig = plt.figure(figsize=(10, 12))  # Adjust the figure size as needed
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])  # Adjust the ratios as per your requirement

    # Magnetogram plot
    ax1 = plt.subplot(gs[0])

    # Magnetic Field Strength and Expansion Factor plot
    # ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2 = ax1.twinx()
    # Make ax2 and ax1 share the same y axis


    # Sunspot Number plot
    ax3 = plt.subplot(gs[1])

    alabel = r"Magnetic Field Strength [Gauss]"
    blabel = "Expansion Factor $F_r$"
    tlabel = r"Latitude [Radians]"

    ax1.set_xlabel(r'Height above Photosphere [r/R$_\odot$-1]')
    ax1.set_ylabel(alabel, color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    # ax1.set_ylim(10**(-3), 10**3)

    # Create a second y-axis for the expansion data
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylim(10**(-4), 10**4)

    ax2.set_ylabel(blabel, color='b', )
    ax2.tick_params(axis='y', labelcolor='b')

    # Load the dat file
    arr = np.loadtxt(filename, skiprows=1).T
    nfluxon = int(arr[0,:].max())

    sunspotplot(ax3)

    # Plot the Curves
    for i in np.arange(0, nfluxon):
        floc = np.where(arr[0,:] == i)[0]
        rr = arr[4, floc]

        expansion = arr[7,floc]
        field_coord = arr[8,floc]
        field_value = arr[9,floc]

        # ax1.plot(rr-1, field_coord, c=plt.cm.Reds(i / nfluxon), label="FLUX bmag" if i == nfluxon-1 else "", alpha=0.8)
        ax1.plot(rr-1, field_value, c=plt.cm.Greens(i / nfluxon), label="FLUX |B|" if i == 3*nfluxon//4 else "", alpha=0.8)
        ax1.plot(rr-1, expansion, c=plt.cm.Blues(i / nfluxon), label=blabel if i == 3*nfluxon//4 else "", alpha=0.8,  zorder=-2000)




    zephyr = load_zephyr()
    ax1.plot(zephyr['rx']-1, zephyr['br'][0], 'r', label='ZEPHYR |B|', zorder = -1001)
    # ax1.plot(zephyr['rx']-1, zephyr['br'].T, cmap=plt.cm.Reds, zorder = -10)
    # Prepare data for LineCollection
    from matplotlib.collections import LineCollection
    n_zep = zephyr['br'].shape[0]
    lines = [np.column_stack([zephyr['rx']-1, zephyr['br'][i]]) for i in range(n_zep)]
    line_segments = LineCollection(lines, cmap='Reds', array=np.arange(n_zep), alpha=0.8, zorder=-1000)

    # Add LineCollection to the plot
    ax1.add_collection(line_segments)


    ax1.legend()
    fig.suptitle(F"Carrington Rotation {args.cr}, with {nfluxon} Open Fields")
    # ax0.set_title(F"Field Strength and Origin at r = {r1:.2f} Rs")
    ax1.set_title(F"Magnetic Field Strength and Expansion Factor")

    # ax0.set_xlabel('Longitude (Radians)')
    # ax0.set_ylabel('Sine latitude')

    # Adjust layout
    ax1.set_xlim(10**(-2.5), 10**2.5)
    ax1.set_ylim(10**(-4), 10**4)
    ax1.grid(True)
    ax2.grid(True)
    fig.set_size_inches((7,7))
    plt.tight_layout()
    # print(" " + bmagname)
    plt.savefig(bmagname)
    if args.show or False:
        plt.show()
    plt.close(fig)
    # print(".", end="", flush=True)
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
            bmagfull, maxlist = multiplot_bmag_all(args, rr, maxlist=maxlist)
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
    assert tried, "Something went wrong with the multiplot_bmag_all function."
    print(f"{args.cr} Done!\n")
    return bmagdir


# Main Code
if __name__ == "__main__":
    # Create the argument parser
    print("\n\tPlotting Bmag_All...", end="")
    from fluxpipe.helpers.pipe_helper import configurations
    configs = configurations()

    parser = argparse.ArgumentParser(description=
                            'This script plots the expansion factor of the given radial_bmag_all.dat')
    parser.add_argument('--cr',     type=int, default=None,    help='Carrington Rotation')
    parser.add_argument('--dat_dir',type=str, default=configs["data_dir"],        help='data directory')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"],      help='select the batch name')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='magnetogram file')
    parser.add_argument('--show',   type=int, default=0)
    parser.add_argument('--file',   type=str, default=None)
    parser.add_argument('--adapt',  type=int, default=configs["adapt"],           help='Use ADAPT magnetograms')
    args = parser.parse_args()


    if len(configs["rotations"]) > 1:
        crlist = configs["rotations"]
        print("crlist", crlist)
        for cr in crlist:
            args.cr = cr
            try:
                bmagdir = run_plots(args)
            except FileNotFoundError:
                print("No file found for CR", cr, "\n")
    else:
        bmagdir = run_plots(args)


    # filename = "expansion_cr{}_f{}.avi".format(args.cr, args.nwant)
    # create_video_from_images(bmagdir, filename)

