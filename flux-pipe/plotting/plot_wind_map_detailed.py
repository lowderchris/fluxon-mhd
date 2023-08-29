"""
Plot the detailed solar wind map
================================

This script is used to plot a detailed solar wind map using various types of plots
including hexagon interpolation, scatter plots, contour plots, and histograms. It takes in
data related to solar wind velocity, magnetic fields, and other parameters to generate these plots.

Attributes:
    --cr: Carrington Rotation (int, default=None)
    --dat_dir: Data directory path (str, default='/Users/cgilbert/vscode/fluxons/fluxon-data')
    --show: (int, default=1)
    --interp: Interpolation method (str, default='linear')
    --nwant: (int, default=0)
    --file: File name to load data from (str, default=None)
    --batch: Batch name (str, default='fluxon_paperfigs_5')

Functions:
    scale_data(): Scale the data between 0 and 1.
    remove_outliers(): Remove outliers from the dataset recursively.
    hex_plot(): Create a hexagon interpolation of the data.
    hist_plot(): Plot a histogram of the velocity data.

"""

# Your code here...


import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
import os.path as path
from scipy.interpolate import griddata

from py_plot_helper import get_ax

from py_pipe_helper import \
    (load_fits_magnetogram, load_magnetogram_params, get_fixed_coords)
from plot_fieldmap import magnet_plot



def scale_data(vel0_clean, vel1_clean, outlier_V0, outlier_V1, scale=15**2, power=1):
    """ Scale the data between 0 and 1, then raise to a power, then scale by a factor

    Parameters
    ----------
    vel0_clean : numpy.ndarray
        The velocity data for the initial state.
    vel1_clean : numpy.ndarray
        The velocity data for the final state.
    outlier_V0 : numpy.ndarray
        The outliers in the initial velocity data.
    outlier_V1 : numpy.ndarray
        The outliers in the final velocity data.
    scale : float, optional
        The scaling factor, by default 15**2.
    power : int, optional
        The power to which the scaled data is raised, by default 1.

    Returns
    -------
    tuple
        Scaled velocity data and scaled outliers for both initial and final states.

    """

    vel0_max = np.nanmax(vel0_clean)
    vel1_max = np.nanmax(vel1_clean)
    vel0_min = np.nanmin(vel0_clean)
    vel1_min = np.nanmin(vel1_clean)

    v0 = scale * ((np.abs(vel0_clean) - vel0_min) / (vel0_max-vel0_min))**power
    v1 = scale * ((np.abs(vel1_clean) - vel1_min) / (vel1_max-vel1_min))**power

    outlier_V0_scaled = scale * ((np.abs(outlier_V0) - vel0_min) / (vel0_max-vel0_min))**power
    outlier_V1_scaled = scale * ((np.abs(outlier_V1) - vel1_min) / (vel1_max-vel1_min))**power

    return v0, v1, outlier_V0_scaled, outlier_V1_scaled

# Remove outliers from the dataset recursively
def remove_outliers(data, ph, th, thresh_low=4, thresh_high=2, n_times= 1):
    """ Remove outliers from the dataset recursively

    Parameters
    ----------
    data : numpy.ndarray
        The data array from which to remove outliers.
    ph : numpy.ndarray
        The phi (azimuthal angle) values corresponding to the data.
    th : numpy.ndarray
        The theta (polar angle) values corresponding to the data.
    thresh_low : int, optional
        The lower threshold for outlier removal, by default 4.
    thresh_high : int, optional
        The upper threshold for outlier removal, by default 2.
    n_times : int, optional
        The number of times to recursively remove outliers, by default 1.

    Returns
    -------
    tuple
        Cleaned data, phi, and theta values, along with the outliers.

    """

    # Get the Good Points
    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    data_good, inds_good = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data)
                        if (mean - thresh_low * std < x < mean + thresh_high * std) and x > 0]))
    ph_good, th_good = ph[inds_good], th[inds_good]
    data_good = np.asarray(data_good)

    # Get the Bad Points
    inds_bad = [i for i in range(len(data)) if i not in inds_good]
    ph_bad, th_bad = ph[inds_bad], th[inds_bad]
    data_bad = data[inds_bad]

    # Return or Recurse
    if n_times <= 1:
        return data_good, ph_good, th_good, data_bad, ph_bad, th_bad
    else:
        data_good_2, ph_good_2, th_good_2, data_bad_2, ph_bad_2, th_bad_2 = \
                remove_outliers(data_good, ph_good, th_good, thresh_low, thresh_high, n_times-1)
        bad_data_all = np.concatenate((data_bad, data_bad_2))
        bad_ph_all   = np.concatenate((  ph_bad,   ph_bad_2))
        bad_th_all   = np.concatenate((  th_bad,   th_bad_2))
        return data_good_2, ph_good_2, th_good_2, bad_data_all, bad_ph_all, bad_th_all


def hex_plot(ph1_clean, th1_clean, vel1_clean, ax=None, nx=20, vmin=400, vmax=800, do_print_top=True, do_hex=True):
    """ Create a hexagon interpolation of the data

    Parameters
    ----------
    ph1_clean : numpy.ndarray
        The phi values of the cleaned data.
    th1_clean : numpy.ndarray
        The theta values of the cleaned data.
    vel1_clean : numpy.ndarray
        The velocity values of the cleaned data.
    ax : matplotlib.pyplot.axis, optional
        An axis to plot on, by default None.
    nx : int, optional
        Number of hexagons in the x direction, by default 20.
    vmin : int, optional
        The minimum value of the colortable, by default 400.
    vmax : int, optional
        The maximum value of the colortable, by default 800.
    do_print_top : bool, optional
        Print with more verbosity, by default True.
    do_hex : bool, optional
        Do the hex plot, else do a contour plot, by default True.

    Returns
    -------
    matplotlib.pyplot.draw_handle
        The output of the plot command.

    """

    if do_print_top:
        if do_hex:
            print("\n\t\tMaking Hexbin Plot...", end="")
        else:
            print("\n\t\tMaking Contour Plot...", end="")
    _, hex_ax = get_ax(ax)
    ## Plot the Interp. data
    # Define the x-boundaries of the domain
    x_min = 0
    x_max = 2*np.pi

    # Wrap the data around the edges of the domain
    ph1 = (ph1_clean)%(2*np.pi)
    th1 = th1_clean
    vel1 = vel1_clean

    ph1_wrapped = np.concatenate((ph1 - x_max, ph1, ph1 + x_max))
    th1_wrapped = np.concatenate((th1, th1, th1))
    vel1_wrapped = np.concatenate((vel1, vel1, vel1))

    # Create a grid for interpolation
    Ny, Nx = load_fits_magnetogram(batch=batch, bo=3, bn=2).shape
    grid_x, grid_y = np.linspace(x_min, x_max, Nx, endpoint=False), np.linspace(-1, 1, Ny)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)

    # Interpolate values on the grid
    points1 = [(ph, th) for ph, th in zip(ph1_wrapped, th1_wrapped)]
    grid_z1 = griddata(points1, vel1_wrapped, (grid_x, grid_y), method=interp, fill_value=0)

    # Plot the interpolated data
    sy, sx = grid_x.shape
    ratio = float(sy)/float(sx)
    gridsize = (int(nx*1.2), int(np.round(nx*ratio)))
    grid_z1_NANed = grid_z1.copy()
    grid_z1_NANed[grid_z1_NANed==0] = np.nan

    z1_use = grid_z1_NANed
    if do_hex:
        hex1 = hex_ax.hexbin(grid_x.flatten(), grid_y.flatten(), C=z1_use.flatten(),
                gridsize=gridsize, cmap='autumn',
                # vmin=np.nanmin(z1_use), vmax=np.nanmax(z1_use))
                vmin=vmin, vmax=vmax)
        print("Done!")
        return hex1
    else:
        # hex1 = hist_ax.imshow(grid_z1, extent=(0, 2*np.pi, -1, 1), zorder=0, alpha=1, cmap="autumn", vmin=vmin, vmax=vmax)
        # contour1 = hex_ax.contourf(grid_z1, zorder=-10, alpha=1, cmap="autumn")
        contour1 = hex_ax.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap="autumn", vmin=vmin, vmax=vmax)
        print("Done!")
        return contour1


def hist_plot(vel1_clean, ax=None, vmin=400, vmax=800, n_bins=20, do_print_top=True):
    """ Plot a histogram of the velocity data

    Parameters
    ----------
    vel1_clean : numpy.ndarray
        The velocity data.
    ax : matplotlib.pyplot.axis, optional
        An axis to plot on, by default None.
    vmin : int, optional
        The minimum value of the colortable, by default 400.
    vmax : int, optional
        The maximum value of the colortable, by default 800.
    n_bins : int, optional
        The resolution of the histogram, by default 20.
    do_print_top : bool, optional
        Print more verbosely, by default True.

    Returns
    -------
    float, float
        Mean and standard deviation of the velocity data.

    """

    if do_print_top:
        print("\n\t\tMaking Histogram Plot...", end='')

    ## The Histogram Plot
    _, hist_ax = get_ax(ax)
    mean1 = np.mean(vel1_clean)
    median1 = np.median(vel1_clean)
    std1  =  np.std(vel1_clean)
    hist_ax.hist(vel1_clean[vel1_clean>=0], bins=n_bins, color='sandybrown')
    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1, label="Mean: {mean1:.0f} km/s")
    hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1, label="Median: {median1:.0f} km/s")
    hist_ax.axvline(mean1+std1, color='k', linestyle=':', linewidth=1, label="Std: {std1:.0f} km/s")
    hist_ax.axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
    hist_ax.legend()
    hist_ax.set_xlabel("Velocity (km/s)")
    hist_ax.set_ylabel("Number of Fluxons")
    hist_ax.set_title(f'CR{CR}, {len(vel1_clean)} Open Fields')

    hist_ax.set_xlim((vmin, vmax))
    if do_print_top:
        print("Success!")
    return mean1, std1
    # ax[3].set_title(f"Solar Wind Velocity Distribution: CR{CR}, Open Fluxons: {len(vel1)}", fontsize=16)

## CODE STARTS HERE
if __name__ == "__main__":
    # create the argument parser
    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxons/fluxon-data', help='data directory')
    parser.add_argument('--show', type=int, default=1)
    parser.add_argument('--interp', type=str, default="linear")
    # parser.add_argument('--nact', type=int, default=0)
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
    vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 3, 1)
    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)
    v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

    ## PLOTTING
    fig, ax = plt.subplots(5)

    mag_ax = ax[0]
    scatter_ax = ax[1]
    contour_ax = ax[2]
    hex_ax = ax[3]
    hist_ax = ax[4]

    all_vmin, all_vmax = 475, 700
    drk=0.25
    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(CR, dat_dir, batch,
        ax=mag_ax, vmin=-500, vmax=500, reduce_amt=reduce, nwant=args.nwant, do_print_top=False)
    hex_n = np.max((n_open//10, 3))

    hex1 = hex_plot(ph1_clean, th1_clean, vel1_clean, ax=hex_ax, nx=hex_n,
                    vmin=all_vmin, vmax=all_vmax)

    scatter_ax.set_facecolor('grey')
    contour_ax.set_facecolor('grey')

    scat1 = scatter_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=6**2, alpha=0.75,
                    marker='s', vmin=all_vmin, vmax=all_vmax, cmap='autumn', edgecolors='none')
    cont1 = contour_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=4**2, alpha=1.0,
                    marker='o', vmin=all_vmin, vmax=all_vmax, cmap='autumn', edgecolors='none')

    scatter_ax.scatter(ph1b, th1b,           color='g', s=3**2, alpha=1, marker='+')
    contour_ax.scatter(ph1b, th1b,           color='g', s=3**2, alpha=1, marker='+')
    scatter_ax.scatter(ph1_clean, th1_clean, c='k', s=2**2, alpha=1., marker='o')


    n_bins = np.linspace(all_vmin, all_vmax, 32)
    mean1, std1 = hist_plot(vel1_clean, ax=hist_ax, vmin=all_vmin, vmax=all_vmax, n_bins=n_bins)


    ## SAVING
    # Set the output file names
    filename = f"png_cr{CR}_f{args.nwant}_ou{n_open}_radial_wind.png"
    main_file =  f'{dat_dir}/batches/{batch}/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    import os
    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))

    for this_ax in [mag_ax, hex_ax, scatter_ax, contour_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.set_ylim((-1.0,1.0))
        this_ax.set_aspect('equal')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline( 1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        this_ax.set_xlim((0, 2*np.pi))

    for jj in [0,1,2]:
        ax[jj].set_xticklabels([])

    # Add an axes for the colorbar
    cbar_ax = fig.add_axes([0.88, 0.2, 0.03, 0.75])

    # Create a colorbar with a custom colormap including the green overlay
    cmap = mpl.colormaps['autumn']

    plotobjs = [scat1, cont1, hex1]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend="both", cmap=cmap, extendfrac=0.1,
                            aspect=15)
        cbar.cmap.set_over('lime')
        cbar.cmap.set_under('darkviolet')
        # cbar.cmap.set_over('lightgreen')
        # cbar.cmap.set_under('lightblue')

    cbar.set_label("Interp. Wind Speed [km/s]", labelpad=-50)

    fig.set_size_inches((6,8))

    plt.subplots_adjust(hspace=0.0,)

    # Save the Figures
    print("\n\t\tSaving figures to disk...", end="")
    # print(main_file)
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace("png", "pdf")

    plt.savefig(outer_file, dpi=200)
    # print("\t\t\tSaving ", shorten_path(outer_pdf))
    plt.savefig(outer_pdf, dpi=200)

    plt.close(fig)
    print("Success!")

    print("\n\t    Done with wind plotting!\n")
    print("\t\t\t```````````````````````````````\n\n\n")
