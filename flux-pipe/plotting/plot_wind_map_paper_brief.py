"""This script plots the solar wind map. It is used to generate plots of the fieldmap and the expansion
"""


import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
import os.path as path
from scipy.interpolate import griddata

from py_plot_helper import get_ax
from py_pipe_helper import load_fits_magnetogram, load_magnetogram_params, get_fixed_coords
from plot_fieldmap import magnet_plot
from plot_wind_map_detailed import hex_plot, hist_plot, scale_data, remove_outliers

## CODE STARTS HERE
if __name__ == "__main__":

    # create the argument parser
    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/fluxons/fluxon-data', help='data directory')
    parser.add_argument('--show', type=int, default=1)
    parser.add_argument('--interp', type=str, default="linear")
    parser.add_argument('--nact', type=int, default=0)
    parser.add_argument('--nwant', type=int, default=0)
    parser.add_argument('--file', type=str, default=None, help='select the file name')
    parser.add_argument('--batch', type=str, default="fluxon_paperfigs_5", help='select the batch name')
    args = parser.parse_args()
    batch = args.batch
    interp = args.interp
    dat_dir = args.dat_dir

    # Load the magnetogram parameters
    (hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(dat_dir)
    CR = args.cr or cr


    # Load the wind file
    dat_file = args.file or f'{dat_dir}/batches/{batch}/cr{CR}/wind/cr{CR}_f{args.nwant}_radial_wind.dat'
    # print("loading file: ", dat_file)
    arr = np.loadtxt(dat_file).T
    try:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
    except ValueError:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
    try:
        nfluxon = arr.shape[1]
    except IndexError as e:
        print("IndexError: ", e)
        print("Wind Calculation Failed, Rerun.")
        exit()


    # Convert coords to correct coords
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)


    print("\n\tPlotting Windmap...", end="\n" if __name__=="__main__" else "")

    # Get the Data
    vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 2, 1)
    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)
    v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

    ## PLOTTING
    fig, ax = plt.subplots(3)

    mag_ax = ax[0]
    hex_ax = ax[1]
    hist_ax = ax[2]

    all_vmin, all_vmax = 450, 850
    drk=0.25

    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(CR, dat_dir, batch,
        ax=mag_ax, vmin=-500, vmax=500, reduce_amt=reduce, nwant=args.nwant, do_print_top=False)
    hex_plot(ph1_clean, th1_clean, vel1_clean, ax=hex_ax, nx=20, vmin=all_vmin, vmax=all_vmax)
    mean1, std1 = hist_plot(vel1_clean, ax=hist_ax, vmin=all_vmin, vmax=all_vmax, n_bins=16)

    sc01 = mag_ax.scatter(ph1_clean, th1_clean, c='k', s=2**2, alpha=0.6, marker='+')
    sc01 = mag_ax.scatter(ph1b, th1b, color=(drk, drk, drk), s=2**2, alpha=0.6, marker='+')


    ## SAVING
    # Set the output file names
    filename = f"cr{CR}_f{args.nwant}_ou{n_open}_radial_wind.png"
    main_file =  f'{dat_dir}/batches/{batch}/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    import os
    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))



    for this_ax in [mag_ax, hex_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.set_ylim((-1.0,1.0))
        this_ax.set_aspect('equal')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline( 1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        this_ax.set_xlim((0, 2*np.pi))

    ax[0].set_xlabel('Longitude (Radians)')
    ax[1].set_xticklabels([])


    fig.set_size_inches((6,6))
    # plt.tight_layout()

    plt.subplots_adjust(
    top=0.981,
    bottom=0.097,
    left=0.113,
    right=0.869,
    hspace=0.312,
    wspace=0.18
    )

    # plt.show()

    # Save the Figures
    print("\n\t\tSaving figures to disk...")
    # print(main_file)
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace(".png", ".pdf")

    plt.savefig(outer_pdf, dpi=200)

    plt.close(fig)
    print("Success!")

    print("\n\t    Done with wind plotting!\n")
    print("\t\t\t```````````````````````````````\n\n\n")
