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
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import shutil
from fluxpipe.helpers.pipe_helper import load_fits_magnetogram, get_fixed_coords
import sunpy.coordinates
import cv2


def plot_bmag_fill(args, r0=0, r1=-1, maxlist=None):
    """
    Plots the expansion factor of the given `radial_bmag_all.dat` file based on the arguments provided.

    Parameters:
        args (argparse.Namespace): An object containing the arguments for plotting.

    Returns:
        True
    """
    configs = configurations(debug=False, args=args)
    batch = args.batch
    ext = "png"
    filename = args.file or f'{args.dat_dir}/batches/{batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    imagename = os.path.basename(filename.replace("all.dat", f"fill_{r0:02d}.{ext}"))
    imagedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(filename))))
    bmagdir = os.path.join(imagedir, "imgs", "bmag", "fill", "sets", f"cr{args.cr}_fill_f{args.nwant}")
    if not os.path.exists(bmagdir):
        os.makedirs(bmagdir)
    bmagname = os.path.join(bmagdir, imagename)


    ### GET TOP AND BOTTOM DATA
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

    low_check = r0
    high_check = r1

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
        maxlist = th0, ph0

    # # Do some data manipulation
    # br0_max = np.nanmax(b0) or 0.25
    # br0_max = np.nanmax(b1) or 0.25
    # # ar0_max = np.nanmax(a0) or 0.25
    # # ar0_max = np.nanmax(a1) or 0.25

    # skew = 5**2
    # power = 1
    # br0 = skew*(6*b0/br0_max)**power
    # br0 = skew*(4*b1/br0_max)**power
    # ar0 = skew*(6*a0/ar0_max)**power
    # ar0 = skew*(4*a1/ar0_max)**power







    # fig, (ax_lowmag, ax_newexpand) = plt.subplots(2)

    RS = 696340000 #meters

    fig2, ax_diff = plt.subplots()

    # Using GridSpec for layout
    fig = plt.figure()  # Adjust the figure size as needed
    gs = gridspec.GridSpec(3, 6, height_ratios=[2.5, 2.5, 1.5])  # Adjust the ratios as per your requirement

    fig.set_size_inches((12,8))


    # Magnetogram plot (top row, spanning two columns)
    ax_lowmag = plt.subplot(gs[0, :3])
    ax_highmag = plt.subplot(gs[0, 3:], sharex=ax_lowmag, sharey=ax_lowmag)
    ax_newexpand = plt.subplot(gs[2, :3])
    ax_cycle = plt.subplot(gs[2, 3:])

    ax_allmag = plt.subplot(gs[1, :2], sharex=ax_lowmag, sharey=ax_lowmag)
    ax_long = plt.subplot(gs[1, 2:4], sharex=ax_lowmag, sharey=ax_newexpand)
    ax_latt = plt.subplot(gs[1, 4:], sharey=ax_newexpand,)

    ax_highmag.tick_params(axis='y', labelleft=False)
    ax_long.tick_params(axis='y', labelleft=False)
    ax_latt.tick_params(axis='y', labelleft=False, labelright=True)



    ### FIRST PLOT ###
    # Plot magnetogram
    magnet = load_fits_magnetogram(configs=configs, cr=args.cr)
    ax_lowmag.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=0.5)
    ax_highmag.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=0.5)
    ax_allmag.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=0.5)

    import warnings
    warnings.filterwarnings('ignore', 'divide by zero encountered in divide')

    #Plot Magnetic Field Strength
    try:
        bb0 = (b0 - np.nanmin(b0)) / (np.nanmax(b0) - np.nanmin(b0))
        if np.isnan(bb0).any():
            bb0[np.isnan(bb0)] = b0[np.isnan(bb0)] / np.nanmax(b0)
    except (ValueError, RuntimeWarning):
        bb0 = b0 / np.nanmax(b0)

    try:
        bb1 = (b1 - np.nanmin(b1)) / (np.nanmax(b1) - np.nanmin(b1))
        if np.isnan(bb1).any():
            bb1[np.isnan(bb1)] = b1[np.isnan(bb1)] / np.nanmax(b1)
    except (ValueError, RuntimeWarning):
        bb1 = b1 / np.nanmax(b1)


    # Plot Latitude of points
    do_lat = True
    if do_lat:
        sc01 = ax_lowmag.scatter(ph0, th0, c=maxlist[0], s = 100, cmap="brg",
                    alpha=bb0, label=r"First Point Latitude", zorder = 99, marker='o', vmin=-1, vmax=1)
        cbar01 = fig.colorbar(sc01, ax=ax_lowmag)
        # cbar01.set_label(f"Lowest Latitude")


    if True:
        sc01 = ax_highmag.scatter(ph1, th1, c=maxlist[0], s = 100, cmap="brg",
                    alpha=bb1, label=r"First Point Latitude", marker='o', vmin=-1, vmax=1)
        cbar01 = fig.colorbar(sc01, ax=ax_highmag)
        cbar01.set_label("Lowest Latitude")




    fig.suptitle(F"Fluxon Magnetic Field Properties for CR {args.cr}, with {nfluxon} Open Fields")
    ax_lowmag.set_title(F"Field Strength and Origin at r = {rad0:.2f} Rs")
    ax_highmag.set_title(F"Field Strength and Origin at r = {rad1:.2f} Rs")
    ax_allmag.set_title(F"Transverse Displacement")
    ax_latt.set_title(F"Lattitude Displacement")
    ax_long.set_title(F"Longitude Displacement")


    # ax_long.set_title(F"Radial Field Characteristics")


    for ax in [ax_lowmag, ax_highmag, ax_allmag]:
        ax.set_xlabel('Longitude [Radians]', labelpad=0)
        ax.set_ylabel('Sine latitude')
        ax.set_aspect('auto')
        ax.set_xlim(-1, 7)
        ax.axvline(0)
        ax.axvline(2*np.pi)
        ax.set_facecolor('grey')



    ax_newexpand.set_title("Fluxon Properties")

    ### SECOND PLOT ###
    fieldlabel = "Magnetic Field Strength"
    expandlabel = r"Expansion Factor"
    latlabel = r"Latitude [Sine Radians]"
    lonlabel = r"Longitude [Radians]"

    ax_latt.set_ylabel(r'Heliocentric Radius [R$_\odot$]')
    ax_latt.set_yscale('log')
    ax_latt.yaxis.set_label_position("right")
    ax_latt.yaxis.set_ticks_position('both')

    ax_newexpand.set_ylabel(r'Heliocentric Radius [R$_\odot$]')
    ax_newexpand.set_yscale('log')
    ax_newexpand.set_xscale('log')

    # Create a second y-axis for the expansion data
    ax_newexpand.set_xscale('log')


    ax_long.set_xlabel(lonlabel)#, color='purple')
    # ax_long.tick_params(axis='x', labelcolor='purple')

    ax_latt.set_xlabel(latlabel)#, color='green')
    # ax_latt.tick_params(axis='x', labelcolor='green')

    radial_heights = []

    # Plot the Curves
    for i in reversed(np.arange(0, nfluxon)):
        # Get data
        floc = np.where(arr[0,:] == i)[0]
        rr = arr[4, floc]

        phi = arr[6, floc]
        ttheta = -1 * arr[5, floc] + np.pi/2
        pphi = np.mod(arr[6, floc], 2*np.pi)


        expansion = arr[7,floc]
        field = arr[8,floc]
        field_2 = arr[9,floc]

        some = 0.1 * nfluxon
        index = (i+some)/(nfluxon+some)
        first_index_lat= (maxlist[0][i]+1)/2
        first_index_long = (maxlist[1][i])/2/np.pi
        col = plt.cm.brg(first_index_lat)



        ppphi = np.unwrap(pphi)
        diff_ppphi = np.abs(np.diff(ppphi))
        diff_rr = np.abs(np.diff(rr))
        dphi_dr = np.abs(diff_ppphi / diff_rr)


        straight_pt = 0.087  # 5 degrees in radians
        #start at one end of the array and iterate until it stops being close to zero, then record that index and plot a line at that index
        for ii in reversed(np.arange(len(dphi_dr))):
            if dphi_dr[ii] > straight_pt:
                break
        ax_diff.scatter(-0.005, rr[ii], color=col, alpha=0.6, zorder=3000, marker="_")
        ax_diff.plot(dphi_dr, rr[:-1], color=col, alpha=0.6, zorder=3000)

        radial_heights.append(rr[ii])

        ### SECOND PLOT ###
        magline, = ax_newexpand.plot(field, rr, c=plt.cm.Reds_r(index), label=fieldlabel if i == 0 else "", alpha=0.7, zorder=2000)
        expline, = ax_newexpand.plot(expansion, rr, c=plt.cm.Blues_r(index), label=expandlabel if i == 0 else "", alpha=0.7, zorder= 1000)

        ### THIRD PLOT ###


        ax_long.plot(ppphi,rr,  color=col, label=latlabel if i == 0 else "", alpha=0.6, zorder=3000-(2*i))
        # ax_long.plot(ppphi,rr,  c=plt.cm.gnuplot(first_index_long), label=latlabel if i == 0 else "", alpha=0.6, zorder=3000)
        ax_long.set_xlim(-1, 2*np.pi+1)

        ### Fourth PLOT ###
        tttheta = np.sin(ttheta)
        ax_latt.plot(tttheta,rr,  color=col, label=lonlabel if i == 0 else "", alpha=0.6, zorder=3000)
        ax_latt.set_xlim(1.1, -1.1)

        ax_allmag.scatter(ppphi[0], tttheta[0], color=col, s=20, linewidth=0, zorder=3000)
        ax_allmag.plot(ppphi, tttheta, color=col, alpha=0.6, zorder=3000)



    # do statistics on the radial heights
    radial_heights = np.array(radial_heights)
    radial_heights = radial_heights[~np.isnan(radial_heights)]
    mean_radial_height = np.mean(radial_heights)
    median_radial_height = np.median(radial_heights)
    std_radial_height = np.std(radial_heights)
    ax_diff.axhline(mean_radial_height, ls="--", c='k', zorder=1000000, label="Mean: {:.2f} R$_\odot$".format(mean_radial_height))
    ax_diff.axhline(median_radial_height, ls=":", c='k', zorder=1000000, label="Median: {:.2f} R$_\odot$".format(median_radial_height))
    ax_diff.axvline(straight_pt, ls=":", c='grey', zorder=1000000, label=f'"Straight": {straight_pt*180/np.pi:.2f} Degrees per R$_\odot$')
    ax_diff.set_xlim(-0.01, 0.1)

    # ax_diff.axhline(mean_radial_height+std_radial_height, ls=":", c='k', zorder=1000000)
    # ax_diff.axhline(mean_radial_height-std_radial_height, ls=":", c='k', zorder=1000000)
    # ax_diff.set_ylim(1, 5)
    ax_diff.set_xlabel(r"$\frac{d\phi}{dr}$ [radians per R$_\odot$]")
    ax_diff.set_ylabel(r"$r$ [R$_\odot$]")
    ax_diff.set_title(f"Radial Expansion, {nfluxon} open fields, Mean: {mean_radial_height:.2f} R$_\odot$")
    # ax_diff.legend()
    # ax_diff.set_title(f"Median Height of Radial Expansion: {median_radial_height:.2f} R$_\odot$, {nfluxon} open fields")
    # ax_diff.set_xlim(0, 1)
    # ax_diff.set_xscale('linear')
    # ax_diff.set_yscale('log')
    # ax_diff.set_xticks([0, 0.5, 1])
    # ax_diff.set_xticklabels([0, 0.5, 1])
    # ax_diff.set_yticks([1, 2, 3, 4, 5])
    # ax_diff.set_yticklabels([1, 2, 3, 4, 5])
    # ax_diff.set_facecolor('grey')
    # ax_diff.set_aspect('auto')
    # ax_diff.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    # ax_diff.grid(True, which='both', axis='both', alpha=0.5, zorder=0)




    r0 = rad0
    r1 = rad1
    ax_newexpand.axhline(r0, ls=":", c='k', zorder=1000000)
    ax_long.axhline(r0, ls=":", c='k', zorder=1000000)
    ax_latt.axhline(r0, ls=":", c='k', zorder=1000000)
    # if r1:
        # ax_newexpand.axvline(r1, ls=":", c='k', zorder=1000)

    # Create handles and labels
    handles, labels = [], []
    for ax in [ax_newexpand]:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        handles.extend(ax_handles)
        labels.extend(ax_labels)

    # Create a single legend
    ax_newexpand.legend(handles, labels, frameon=False, loc='center left')
    ax_newexpand.set_xlim((1* 10**(-6), 5* 10**2))




    ### Fourth PLOT ###
    # Plot the Sunspot Number
    carrington = np.loadtxt("/Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/SN_m_tot_V2.0.txt").T
    ## https://sidc.be/SILSO/datafiles#total ##
    import sunpy.coordinates
    date = carrington[2]
    sunspots = carrington[3]
    this_date = sunpy.coordinates.sun.carrington_rotation_time(args.cr)
    # CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))

    # fig, ax = plt.subplots()
    ax_cycle.plot(date, sunspots, label="Sunspots")
    ax_cycle.axvline(this_date.decimalyear, ls=":", c='k', zorder=1000000)
    ax_cycle.set_xlabel("Year")
    ax_cycle.set_ylabel("Sunspots")
    ax_cycle.yaxis.set_label_position("right")
    ax_cycle.yaxis.set_ticks_position('right')

    ax_cycle.set_title("Solar Cycle")
    ax_cycle.set_xlim(2005, 2025)
    ax_cycle.set_ylim(0, 200)
    # fig.close()
    # plt.show()

    # ax_lowmag.set_xlim(0, 2*np.pi)





    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(
        top=0.92,
        bottom=0.063,
        left=0.063,
        right=0.95,
        hspace=0.365,
        wspace=0.1
    )


    fig.savefig(bmagname)

    dr = os.path.dirname(os.path.dirname(os.path.dirname(bmagname))) + "/radial"
    if not os.path.exists(dr):
        os.makedirs(dr)
    pathh = os.path.join(dr, os.path.basename(bmagname).replace(".png", "_radial.png"))
    print(pathh)
    fig2.savefig(pathh)

    if args.show or False:
        plt.show()
    plt.close(fig)
    plt.close(fig2)
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

    video_folder = os.path.join(os.path.dirname(image_folder), "..", "expansion_videos")
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
    bmagdir=None
    for rr in np.arange(times+1):
        tried = True
        try:
            bmagfull, maxlist = plot_bmag_fill(args, rr, maxlist=maxlist)
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
        assert tried, "Something went wrong with the plot_bmag_fill function."
    print("Done!\n")
    return bmagdir


# Main Code
if __name__ == "__main__":
    # Create the argument parser
    print("\n\tPlotting Bmag_Fill...", end="")
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
    create_video_from_images(bmagdir, filename)

