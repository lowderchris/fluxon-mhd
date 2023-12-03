"""
Plot the detailed solar wind map
================================

This script is used to plot a detailed solar wind map using various types of plots
including hexagon interpolation, scatter plots, contour plots, and histograms. It takes in
data related to solar wind velocity, magnetic fields, and other parameters to generate these plots.

Attributes:
    --cr: Carrington Rotation (int, default=None)
    --dat_dir: Data directory path (str, default is defined in config.ini)
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


import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
import os

import os.path as path
from scipy.interpolate import griddata
from scipy.stats import norm
from scipy.optimize import curve_fit
from fluxpipe.helpers.pipe_helper import (configurations, load_fits_magnetogram,
                         load_magnetogram_params, get_fixed_coords, get_ax)

from fluxpipe.plotting.plot_fieldmap import magnet_plot


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

    if n_times == 0:
        return data, ph, th, np.array([]), np.array([]), np.array([])
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


def hex_plot(ph1_clean, th1_clean, vel1_clean, ax=None, nx=20, vmin=400, vmax=800, do_print_top=True, do_hex=False, configs=None, CR=None, cmap='autumn'):
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
    Ny, Nx = load_fits_magnetogram(configs=configs, bo=3, bn=2).shape
    grid_x, grid_y = np.linspace(x_min, x_max, Nx, endpoint=False), np.linspace(-1, 1, Ny)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)

    # Interpolate values on the grid
    points1 = [(ph, th) for ph, th in zip(ph1_wrapped, th1_wrapped)]
    grid_z1 = griddata(points1, vel1_wrapped, (grid_x, grid_y), method=args.interp, fill_value=0)

    # Plot the interpolated data
    sy, sx = grid_x.shape
    ratio = float(sy)/float(sx)
    gridsize = (int(nx*1.2), int(np.round(nx*ratio)))
    grid_z1_NANed = grid_z1.copy()
    grid_z1_NANed[grid_z1_NANed==0] = np.nan

    z1_use = grid_z1_NANed
    if do_hex:
        hex1 = hex_ax.hexbin(grid_x.flatten(), grid_y.flatten(), C=z1_use.flatten(),
                gridsize=gridsize, cmap=cmap,
                # vmin=np.nanmin(z1_use), vmax=np.nanmax(z1_use))
                vmin=vmin, vmax=vmax)
        print("Success!")
        return hex1
    else:
        hex1 = hex_ax.imshow(grid_z1, extent=(0, 2*np.pi, -1, 1), zorder=0, alpha=1, cmap=cmap, vmin=vmin, vmax=vmax)
        # contour1 = hex_ax.contourf(grid_z1, zorder=-10, alpha=1, cmap="autumn")
        # contour1 = hex_ax.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap=cmap, vmin=0.5*vmin, vmax=2*vmax)
        print("Success!")
        return hex1


def hist_plot_orig(vel1_clean, ax=None, vmin=400, vmax=800, n_bins=20, do_print_top=True, CR="<unset>"):
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
    fig, hist_ax = get_ax(ax)
    mean1 = np.mean(vel1_clean)
    median1 = np.median(vel1_clean)
    std1  =  np.std(vel1_clean)
    vel_positive = vel1_clean[vel1_clean>=0]
    hist_ax.hist(vel_positive, bins=n_bins, color='sandybrown')
    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1, label=f"Mean: {mean1:.0f} km/s")
    hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1, label=f"Median: {median1:.0f} km/s")
    hist_ax.axvline(mean1+std1, color='k', linestyle=':', linewidth=1, label=f"Std: {std1:.0f} km/s")
    hist_ax.axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
    hist_ax.legend()
    hist_ax.set_xlabel("Velocity (km/s)")
    hist_ax.set_ylabel("Number of Fluxons")
    fig.suptitle(f'CR{CR}, {len(vel1_clean)} Open Fields')

    hist_ax.set_xlim((vmin, vmax))
    if do_print_top:
        print("Success!")

    return mean1, std1


def hist_plot(vel1_clean, ax=None, vmin=400, vmax=800, n_bins=20, do_print_top=True,
              CR="<unset>", cmap='autumn', do_gaussian = True):
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
    CR : str, optional
        Carrington Rotation information, by default "<unset>".

    Returns
    -------
    float, float
        Mean and standard deviation of the velocity data.

    """

    if do_print_top:
        print("\n\t\tMaking Histogram Plot...", end='')

    ## The Histogram Plot
    fig, hist_ax = get_ax(ax)
    mean1 = np.mean(vel1_clean)
    median1 = np.median(vel1_clean)
    std1 = np.std(vel1_clean)
    vel_positive = vel1_clean[vel1_clean >= 0]

    # Calculate the histogram
    hist, bin_edges = np.histogram(vel_positive, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot the histogram
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    for i in range(len(bin_edges) - 1):
        color = cmap(norm(bin_centers[i]))
        hist_ax.bar(
            bin_edges[i], hist[i], width=bin_edges[i + 1] - bin_edges[i], color=color
        )

    if do_gaussian:
        # Define the Gaussian (normal) distribution function
        def gaussian(x, amplitude, mean, stddev, b):
            return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2) + b

        # Fit the Gaussian curve to the histogram data
        params, _ = curve_fit(gaussian, bin_centers, hist, p0=[1, np.median(vel_positive), np.std(vel_positive),0])

        # Extract the fitted parameters
        amplitude, mean, stddev, b = params

        # Generate the fitted curve
        fitted_curve = gaussian(bin_centers, *params)

        plt.plot(bin_centers, fitted_curve, 'b-', label='Fitted Gaussian Curve')

        mean1 = mean
        std1 = np.abs(stddev)

        mean_label = f"GaussMean: {mean1:.0f} km/s"
        median_label = f"StatMedian: {median1:.0f} km/s"
        std_label = f"GaussStd: {std1:.0f} km/s"
        # b_label = f"GaussB: {b:.0f}"
        # hist_ax.axhline(b, color='gray', linestyle='-', linewidth=1, label=b_label)
    else:
        mean_label = f"StatMean: {mean1:.0f} km/s"
        median_label = f"StatMedian: {median1:.0f} km/s"
        std_label = f"StatStd: {std1:.0f} km/s"


    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1, label=mean_label)
    hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1, label=median_label)
    hist_ax.axvline(mean1 + std1, color='k', linestyle=':', linewidth=1, label=std_label)
    hist_ax.axvline(mean1 - std1, color='k', linestyle=':', linewidth=1)

    hist_ax.legend(fontsize="small")
    hist_ax.set_xlabel("Velocity (km/s)")
    hist_ax.set_ylabel("Number of Fluxons")
    fig.suptitle(f'CR{CR}, {len(vel1_clean)} Open Fields')

    hist_ax.set_xlim((vmin-50, vmax+100))

    if do_print_top:
        print("Success!")

    return mean1, std1

def plot_wind_map_detailed_orig(configs):
    """ Plot the detailed solar wind map

    Parameters
    ----------
    args : argparse.Namespace
        The arguments to the script.

    Returns
    -------
    bool
        True if the plot was successful, else False.

    """

    print("\n\tPlotting Windmap...", end="\n" if __name__=="__main__" else "")
    configs = configs or configurations()

    # Set the arguments
    batch =     configs.get("batch_name")
    dat_dir =   configs.get("data_dir")
    nwant =     configs.get("nwant")
    CR =        configs.get("cr")
    dat_file =  configs.get("file", f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/cr{CR}_f{nwant}_radial_wind.dat')

    # Load the wind file
    arr = np.loadtxt(dat_file).T

    try:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
    except ValueError:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
    try:
        nfluxon = arr.shape[1]
    except IndexError as e:
        print("IndexError: ", e)
        print("Not enough open field lines. Increase Fluxon count.")
        exit()

    # Convert coords to correct coords
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)


    # Get the Data
    # vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 3, 1)
    # v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)
    # vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 100, 100, 1)
    # vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 100, 100, 1)

    ## PLOTTING
    fig, ax = plt.subplots(5)

    mag_ax = ax[0]
    square_ax = ax[1]
    dot_ax = ax[2]
    hex_ax = ax[3]
    hist_ax = ax[4]

    all_vmin, all_vmax = 450, 700
    drk=0.25

    cmap = mpl.colormaps["autumn"]
    cmap.set_over('lime')
    cmap.set_under('darkviolet')


    # if configs["adapt"]:
    #     reduce_amt = 'A'
    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(configs=configs,
        ax=mag_ax, vmin=-500, vmax=500, do_print_top=False, legend=True)
    mag_ax.set_title("R=Rs Magnetogram with Footpoints")
    mag_ax.set_xlabel("Longitude [rad]", labelpad=-3)


    hex_n = np.max((n_open//10, 3))

    hex1 = hex_plot(ph1_clean, th1_clean, vel1_clean, ax=hex_ax, nx=hex_n,
                    vmin=all_vmin, vmax=all_vmax, configs=configs, cmap=cmap)

    hex_ax.set_title("R=21Rs Wind Speed")
    square_ax.set_title("R=21Rs Wind Speed")
    dot_ax.set_title("R=21Rs Wind Speed")

    square_ax.set_facecolor('grey')
    dot_ax.set_facecolor('grey')

    scat1 = square_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=6**2, alpha=0.75,
                    marker='s', vmin=all_vmin, vmax=all_vmax, cmap=cmap, edgecolors='none')
    cont1 = dot_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=4**2, alpha=1.0,
                    marker='o', vmin=all_vmin, vmax=all_vmax, cmap=cmap, edgecolors='none')

    square_ax.scatter(ph1b, th1b,           color='w', s=3**2, alpha=1, marker='+')
    dot_ax.scatter(ph1b, th1b,           color='w', s=3**2, alpha=1, marker='+')
    # square_ax.scatter(ph1_clean, th1_clean, c='k', s=2**2, alpha=1., marker='o')


    n_bins = np.linspace(all_vmin//2, 2*all_vmax, 64)
    mean1, std1 = hist_plot(vel1_clean, ax=hist_ax, vmin=all_vmin, vmax=all_vmax, n_bins=n_bins, CR=CR, cmap=cmap)
    hist_ax.set_title("Histogram of Wind Speeds at 21Rs")






    ## SAVING
    # Set the output file names
    filename = f"png_cr{CR}_f{nwant}_ou{n_open}_radial_wind.png"
    main_file =  f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))

    for this_ax in [mag_ax, hex_ax, square_ax, dot_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline( 1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        this_ax.set_aspect('equal')
        this_ax.set_ylim((-1.0,1.0))
        this_ax.set_xlim((0, 2*np.pi))

    for jj in [0,1,2]:
        ax[jj].set_xticklabels([])

    # Add an axes for the colorbar
    cbar_ax = fig.add_axes([0.9, -0.085, 0.03, 1.0])

    plotobjs = [scat1, cont1, hex1]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend="both", cmap=cmap, extendfrac=0.1, aspect=15)


    cbar.set_label("Interpolated Wind Speed [km/s]", labelpad=-50)

    fig.set_size_inches((6,10))

    plt.subplots_adjust(top=0.940,
                        bottom=0.065,
                        left=0.12,
                        right=0.865,
                        hspace=0.340,
                        wspace=0.155)

    # Save the Figures
    print("\n\t\tSaving figures to disk...", end="")
    # print(main_file)
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace("png", "pdf")

    plt.show()
    plt.savefig(outer_file, dpi=200)
    # print("\t\t\tSaving ", shorten_path(outer_pdf))
    # plt.savefig(outer_pdf, dpi=200)

    plt.close(fig)
    print("Success!")

    print("\n\t    Done with wind plotting!\n")
    print("\t\t\t```````````````````````````````\n\n\n")


def plot_wind_map_detailed(configs):
    """
    Plot the detailed solar wind map.

    Parameters
    ----------
    configs : dict
        Dictionary containing configuration settings.

    Returns
    -------
    bool
        True if the plot was successful, else False.
    """
    print("\n\tPlotting Windmap...", end="\n" if __name__ == "__main__" else "")
    configs = configs or configurations()

    # Extract configuration settings
    batch = configs.get("batch_name")
    dat_dir = configs.get("data_dir")
    nwant = configs.get("nwant")
    CR = configs.get("cr")
    dat_file = configs.get("file", f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/cr{CR}_f{nwant}_radial_wind.dat')

    all_vmin, all_vmax = configs.get("all_vmin", 450), configs.get("all_vmax", 700)

    # Load the wind file
    arr = np.loadtxt(dat_file).T

    try:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
    except ValueError:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
    try:
        nfluxon = arr.shape[1]
    except IndexError as e:
        print("IndexError: ", e)
        print("Not enough open field lines. Increase Fluxon count.")
        exit()

    # Convert coordinates to the correct format
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)

    # Clean and process data (remove outliers, scale, etc.)
    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)

    n_total = len(vel1)
    n_clean = len(vel1_clean)

    n_outliers = len(vel1[vel1<all_vmin]) + len(vel1[vel1>all_vmax])
    percent_outliers = 100 * n_outliers / n_total


    # Create subplots
    fig, ax = plt.subplots(5)

    mag_ax, square_ax, dot_ax, hex_ax, hist_ax = ax

    cmap = mpl.colormaps["autumn"]
    cmap.set_over('lime')
    cmap.set_under('darkviolet')

    # Plot R=Rs Magnetogram with Footpoints
    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(
        configs=configs,
        ax=mag_ax,
        vmin=-500,
        vmax=500,
        do_print_top=False,
        legend=True
    )


    # Interpolated plot of wind speed [Hexagons]
    hex_n = np.max((n_open // 10, 3))
    hex1 = hex_plot(
        ph1_clean,
        th1_clean,
        vel1_clean,
        do_hex=False,
        ax=hex_ax,
        nx=hex_n,
        vmin=all_vmin,
        vmax=all_vmax,
        configs=configs,
        cmap=cmap
    )

    # Scatter plot of wind speed [Squares]
    scat1 = square_ax.scatter(
        ph1_clean,
        th1_clean,
        c=vel1_clean,
        s=6**2,
        alpha=0.75,
        marker='s',
        vmin=all_vmin,
        vmax=all_vmax,
        cmap=cmap,
        edgecolors='none'
    )

    # Scatter plot of wind speed [Circles]
    cont1 = dot_ax.scatter(
        ph1_clean,
        th1_clean,
        c=vel1_clean,
        s=4**2,
        alpha=1.0,
        marker='o',
        vmin=all_vmin,
        vmax=all_vmax,
        cmap=cmap,
        edgecolors='none'
    )

    # Plot the Outliers
    square_ax.scatter   (ph1b, th1b, color='w', s=3**2, alpha=1, marker='+')
    dot_ax.scatter      (ph1b, th1b, color='w', s=3**2, alpha=1, marker='+')


    # Create histogram of wind speeds
    n_bins = np.linspace(all_vmin - 200, all_vmax + 300, 48)
    mean1, std1 = hist_plot(
        vel1_clean,
        ax=hist_ax,
        vmin=all_vmin,
        vmax=all_vmax,
        n_bins=n_bins,
        CR=CR,
        cmap=cmap
    )


    ## Plot Formatting #####
    # Set Titles
    hist_ax.set_title("Histogram of Wind Speeds at r=21Rs")
    mag_ax.set_title("Magnetogram at r=Rs, with Footpoints")
    hex_ax.set_title("Wind Speed at R=21Rs")
    square_ax.set_title(f"Wind Speed at R=21Rs, with {percent_outliers:0.1f}% Outliers")
    dot_ax.set_title("Wind Speed at R=21Rs")

    mag_ax.set_xlabel("Longitude [rad]", labelpad=-3)

    square_ax.set_facecolor('grey')
    dot_ax.set_facecolor('grey')

    for this_ax in [mag_ax, hex_ax, square_ax, dot_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline(1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2 * np.pi, c='lightgrey', zorder=-10)
        this_ax.set_aspect('equal')
        this_ax.set_ylim((-1.0, 1.0))
        this_ax.set_xlim((0, 2 * np.pi))

    for jj in [0, 1, 2]:
        ax[jj].set_xticklabels([])

    # Add a colorbar
    cbar_ax = fig.add_axes([0.9, -0.085, 0.03, 1.0])

    plotobjs = [scat1, cont1, hex1]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend="both", cmap=cmap, extendfrac=0.1, aspect=15)

    cbar.set_label("Interpolated Wind Speed [km/s]", labelpad=-50)

    fig.set_size_inches((6, 10))

    plt.subplots_adjust(
        top=0.940,
        bottom=0.065,
        left=0.12,
        right=0.865,
        hspace=0.340,
        wspace=0)


    # Set the output file names
    filename = f"png_cr{CR}_f{nwant}_ou{n_open}_radial_wind.png"
    main_file = f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))


    # Save the Figures
    print("\n\t\tSaving figures to disk...", end="")
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace("png", "pdf")

    # plt.show()
    plt.savefig(outer_file, dpi=200)
    plt.close(fig)
    print("Success!")

    print("\nDone with wind plotting!\n")
    print("```````````````````````````````\n\n\n")


########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    # Create the argument parser
    configs = configurations()

    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr',     type=int, default=configs["rotations"][0], help='Carrington Rotation')
    parser.add_argument('--interp', type=str, default="linear")
    parser.add_argument('--show',   type=int, default=configs["verbose"])
    parser.add_argument('--file',   type=str, default=None, help='select the file name')
    parser.add_argument('--dat_dir',type=str, default=configs['data_dir'], help='data directory')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='number of fluxons')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"], help='select the batch name')
    parser.add_argument('--adapt',  type=int, default=configs["adapt"], help='Use ADAPT magnetograms')
    args = parser.parse_args()
    configs = configurations(args=args)


    plot_wind_map_detailed(configs)
