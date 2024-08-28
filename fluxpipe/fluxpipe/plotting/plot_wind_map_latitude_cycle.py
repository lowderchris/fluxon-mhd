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
from fluxpipe.science.pfss_funcs import get_fluxon_locations

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
def remove_outliers(data, ph, th, fr, thresh_low=4, thresh_high=2, n_times= 1):
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
        return data, ph, th, fr, np.array([]), np.array([]), np.array([]), np.array([])
    # Get the Good Points
    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    data_good, inds_good = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data)
                        if (mean - thresh_low * std < x < mean + thresh_high * std) and x > 0]))
    ph_good, th_good, fr_good = ph[inds_good], th[inds_good], fr[inds_good]
    data_good = np.asarray(data_good)

    # Get the Bad Points
    inds_bad = [i for i in range(len(data)) if i not in inds_good]
    ph_bad, th_bad = ph[inds_bad], th[inds_bad]
    data_bad = data[inds_bad]

    # Return or Recurse
    if n_times <= 1:
        return data_good, ph_good, th_good, fr_good, data_bad, ph_bad, th_bad, inds_good
    else:
        data_good_2, ph_good_2, th_good_2, data_bad_2, ph_bad_2, th_bad_2 = \
                remove_outliers(data_good, ph_good, th_good, thresh_low, thresh_high, n_times-1)
        bad_data_all = np.concatenate((data_bad, data_bad_2))
        bad_ph_all   = np.concatenate((  ph_bad,   ph_bad_2))
        bad_th_all   = np.concatenate((  th_bad,   th_bad_2))
        return data_good_2, ph_good_2, th_good_2, fr_good, bad_data_all, bad_ph_all, bad_th_all, inds_good


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
    _, high_lat_ax = get_ax(ax)
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
        hex1 = high_lat_ax.hexbin(grid_x.flatten(), grid_y.flatten(), C=z1_use.flatten(),
                gridsize=gridsize, cmap=cmap,
                # vmin=np.nanmin(z1_use), vmax=np.nanmax(z1_use))
                vmin=vmin, vmax=vmax)
        print("Success!")
        return hex1
    else:
        hex1 = high_lat_ax.imshow(grid_z1, extent=(0, 2*np.pi, -1, 1), zorder=0, alpha=1, cmap=cmap, vmin=vmin, vmax=vmax, origin="lower")
        # contour1 = high_lat_ax.contourf(grid_z1, zorder=-10, alpha=1, cmap="autumn")
        # contour1 = high_lat_ax.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap=cmap, vmin=0.5*vmin, vmax=2*vmax)
        print("Success!")
        return hex1, grid_z1


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
        try:
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
            # median_label = f"StatMedian: {median1:.0f} km/s"
            std_label = f"GaussStd: {std1:.0f} km/s"
            # b_label = f"GaussB: {b:.0f}"
            # hist_ax.axhline(b, color='gray', linestyle='-', linewidth=1, label=b_label)
        except:
            mean_label = f"StatMean: {mean1:.0f} km/s"
            # median_label = f"StatMedian: {median1:.0f} km/s"
            std_label = f"StatStd: {std1:.0f} km/s"
    else:
        mean_label = f"StatMean: {mean1:.0f} km/s"
        # median_label = f"StatMedian: {median1:.0f} km/s"
        std_label = f"StatStd: {std1:.0f} km/s"


    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1, label=mean_label)
    # hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1, label=median_label)
    hist_ax.axvline(mean1 + std1, color='k', linestyle=':', linewidth=1, label=std_label)
    hist_ax.axvline(mean1 - std1, color='k', linestyle=':', linewidth=1)

    hist_ax.legend(fontsize="small", loc = "center left")
    hist_ax.set_xlabel("Wind Speed [km/s]")
    hist_ax.set_ylabel("Number of Fluxons")
    fig.suptitle(f'CR{CR}, {len(vel1_clean)} Open Fields, {configs.get("flow_method").title()} wind')

    hist_ax.set_xlim((vmin-25, vmax+25))

    if do_print_top:
        print("Success!")

    return mean1, std1

def plot_closed_footpoints(ax, closed_f):
    """
    Scatter plots closed magnetic footpoints on the given pyplot axis.

    Parameters:
    - ax: Matplotlib axis object where footpoints will be plotted.
    - closed_f: Filepath to the data for closed magnetic footpoints.
    """

    all_file    = closed_f.replace("closed_", "")
    f_lat, f_lon, f_sgn, _fnum = get_fluxon_locations(all_file, None, configs=None)

    # Filter positive and negative cases
    positive_indices = [i for i, s in enumerate(f_sgn) if s > 0]
    negative_indices = [i for i, s in enumerate(f_sgn) if s <= 0]

    # Data for positive and negative cases
    f_lon_positive = [f_lon[i] for i in positive_indices]
    f_lat_positive = [f_lat[i] for i in positive_indices]
    f_lon_negative = [f_lon[i] for i in negative_indices]
    f_lat_negative = [f_lat[i] for i in negative_indices]

    # Plot positive cases with labels
    ax.scatter(f_lon_positive, f_lat_positive, s=3**2, c='lightcyan', alpha=0.4, label='Positive', edgecolors='none')

    # Plot negative cases with labels
    ax.scatter(f_lon_negative, f_lat_negative, s=3**2, c='bisque', alpha=0.4, label='Negative', edgecolors='none')

def plot_wind_map_latitude(configs):
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
    configs = configs or configurations()

    # Extract configuration settings
    batch = configs.get("batch_name")
    dat_dir = configs.get("data_dir")
    nwant = configs.get("nwant")
    reduce = configs.get("mag_reduce")
    CR = configs.get("cr")
    dat_file = configs.get("file", f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/cr{CR}_f{nwant}_radial_wind.dat')
    all_vmin, all_vmax = configs.get("all_vmin", 450), configs.get("all_vmax", 700)
    print(f"\n\tPlotting Windmap {CR} Lat Cycle...", end="\n" if __name__ == "__main__" else "")

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


    n_open = len(phi0)

    # Convert coordinates to the correct format
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)

    ph0 += np.pi
    ph1 += np.pi
    ph0 = np.mod(ph0, 2*np.pi)
    ph1 = np.mod(ph1, 2*np.pi)
    # ph1 = np.unwrap(ph1)

    # Clean and process data (remove outliers, scale, etc.)
    vel1_clean, ph1_clean, th1_clean, fr1_clean, v1b, ph1b, th1b, inds_good = remove_outliers(vel1, ph1, th1, fr1, 3, 2, 1)

    n_total = len(vel1)
    n_clean = len(vel1_clean)

    n_outliers = len(vel1[vel1<all_vmin]) + len(vel1[vel1>all_vmax])
    percent_outliers = 100 * n_outliers / n_total


    pole_file = f'{dat_dir}/batches/{batch}/data/cr{CR}/floc/floc_open_cr{CR}_r{reduce}_f{nwant}_hmi.dat'
    arr2 = np.loadtxt(pole_file).T
    ffid, pol, _, _, _ = arr2
    # Get unique values of fid
    unique_fid = np.unique(ffid)
    # Initialize an empty list to store the result
    polar = []
    # Iterate over unique values of fid
    for unique_value in unique_fid:
        index = np.where(ffid == unique_value)[0][0]
        polar.append(pol[index])
    polarity = np.array(polar)




    import matplotlib.gridspec as gridspec
    # Create a 6x1 grid with custom height ratios
    fig = plt.figure(figsize=(6, 8))
    height_ratios = [ 1, 1, 0.4, 1, 1, 1]  # Adjust the height ratios as needed
    gs = gridspec.GridSpec(6, 1, height_ratios=height_ratios)
    # Create subplots using gridspec
    carr_ax =       plt.subplot(gs[2])
    low_lat_ax =    plt.subplot(gs[0])
    high_lat_ax =   plt.subplot(gs[1])
    square_wind_ax =plt.subplot(gs[3])
    dot_wind_ax =   plt.subplot(gs[4])
    hist_ax =       plt.subplot(gs[5])






    # Plot the Sunspot Number
    carrington = np.loadtxt("/Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/SN_m_tot_V2.0.txt").T
    ## https://sidc.be/SILSO/datafiles#total ##
    import sunpy.coordinates
    date = carrington[2]
    sunspots = carrington[3]
    this_date = sunpy.coordinates.sun.carrington_rotation_time(args.cr)
    # CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))

    # fig, ax = plt.subplots()
    carr_ax.plot(date, sunspots, label="Sunspots")
    carr_ax.axvline(this_date.decimalyear, ls=":", c='k', zorder=1000000)
    # carr_ax.set_xlabel("Year")
    carr_ax.set_ylabel("Sunspots")
    carr_ax.set_title("Solar Cycle Phase")
    carr_ax.set_xlim(2005, 2025)
    # set the major tick formatter to display integers
    from matplotlib.ticker import MaxNLocator
    carr_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    carr_ax.set_ylim(0, 200)









    cmap = mpl.colormaps["autumn"]
    # cmap.set_over('lime')
    # cmap.set_under('darkviolet')

    # # Plot R=Rs Magnetogram with Footpoints


    # Plot magnetogram
    magnet = load_fits_magnetogram(configs=configs, cr=args.cr)
    low_lat_ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=1)

    high_lat_ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
            extent=(0,2*np.pi,-1,1), aspect='auto', vmin=-500, vmax=500, alpha=1)

    the_cmap = "Spectral"

    sc01 = low_lat_ax.scatter(ph0, th0, c=th0, s = 60, cmap=the_cmap, edgecolors='k', linewidths=0.25,
                alpha=0.75, zorder = 100, marker='o', vmin=-1, vmax=1)


    # Add a colorbar
    cbar_ax_lat1 = fig.add_axes([0.9, 0.815, 0.03, 0.123])
    cbar01 = plt.colorbar(sc01, cax=cbar_ax_lat1, cmap=the_cmap, aspect=15)
    cbar01.set_label(f"Base Latitude", labelpad=-58)



    sc02 = high_lat_ax.scatter(ph1, th1, c=th0, s = 60, cmap=the_cmap, zorder=100, edgecolors='k', linewidths=0.25,
                alpha=0.75, label="First Point Latitude", marker='o', vmin=-1, vmax=1)

    # Add a colorbar
    cbar_ax_lat2 = fig.add_axes([0.9, 0.651, 0.03, 0.123])
    cbar02 = plt.colorbar(sc02, cax=cbar_ax_lat2, cmap=the_cmap, aspect=15)
    cbar02.set_label(f"Base Latitude", labelpad=-58)



    nwant = configs.get("nwant", None)
    reduce_amt = configs.get("mag_reduce", None)
    if configs['adapt']:
        inst = "adapt"
        reduce_amt = "f" + str(configs.get("adapt_select"))
    else:
        inst = "hmi"

    floc_path = f"{dat_dir}/batches/{batch}/data/cr{CR}/floc/"
    closed_file = f"{floc_path}floc_closed_cr{CR}_r{reduce_amt}_f{nwant}_{inst}.dat"
    plot_closed_footpoints(low_lat_ax, closed_file)





    # Interpolated plot of wind speed [Hexagons]
    hex_n = np.max((n_open // 40, 3))
    hex1, grid_z1 = hex_plot(
        ph1_clean,
        th1_clean,
        vel1_clean,
        do_hex=False,
        ax=dot_wind_ax,
        nx=hex_n,
        vmin=all_vmin,
        vmax=all_vmax,
        configs=configs,
        cmap=cmap
    )

    clr = np.log(2*fr1/fr0)**1/2
    clr = np.log(fr1/fr0)
    sc11 = square_wind_ax.scatter(ph1, th1, c=clr, s = (4*fr1/fr0)**1/2, cmap="YlGnBu", alpha=0.75, label=r"A(21.5Rs)", marker='s', vmin=-0.3, vmax=0)

    # Add a colorbar
    cbar_ax_lat3 = fig.add_axes([0.91, 0.395, 0.03, 0.137])
    cbar03 = plt.colorbar(sc11, cax=cbar_ax_lat3, aspect=15, extend="max", extendfrac=0.1)
    cbar03.set_label(f"log(Expansion Ratio)", labelpad=-49)
    # # Scatter plot of wind speed [Squares]
    # scat1 = square_wind_ax.scatter(
    #     ph1_clean,
    #     th1_clean,
    #     c=vel1_clean,
    #     s=6**2,
    #     alpha=0.75,
    #     marker='s',
    #     vmin=all_vmin,
    #     vmax=all_vmax,
    #     cmap=cmap,
    #     edgecolors='none'
    # )
#
    # Scatter plot of wind speed [Circles]
    cont1 = dot_wind_ax.scatter(
        ph1_clean,
        th1_clean,
        c=vel1_clean,
        s=4**2,
        lw = 0.1,
        alpha=0.6,
        marker='o',
        vmin=all_vmin,
        vmax=all_vmax,
        cmap=cmap,
        edgecolors='gray',
        # edgewidth=0.5
    )
    # cont1 = None
    # Plot the Outliers
    # square_wind_ax.scatter   (ph1b, th1b, color='lightgrey', s=3**2, alpha=1, marker='+')
    # dot_wind_ax.scatter      (ph1b, th1b, color='lightgrey', s=3**2, alpha=1, marker='+')


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
    hist_ax.set_title("Histogram of Wind Speeds at r=21Rs", pad=-3)
    low_lat_ax.set_title("Fluxon Locations at R=Rs", pad=-3)
    high_lat_ax.set_title("Fluxon Locations at R=21.5Rs", pad=-3)
    dot_wind_ax.set_title(f"Wind Speed at R=21.5Rs, with {percent_outliers:0.1f}% Outliers", pad=-3)
    square_wind_ax.set_title("Expansion Factor Ratio log( fr(21.5Rs)/fr(Rs) )", pad=-3)

    low_lat_ax.set_xlabel("Longitude [rad]", labelpad=-3)

    square_wind_ax.set_facecolor('grey')
    dot_wind_ax.set_facecolor('grey')
    ax_list = [low_lat_ax, high_lat_ax, square_wind_ax, dot_wind_ax]

    for this_ax in ax_list:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline(1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2 * np.pi, c='lightgrey', zorder=-10)
        this_ax.set_aspect('equal')
        this_ax.set_ylim((-1.0, 1.0))
        this_ax.set_xlim((0, 2 * np.pi))
        this_ax.grid(True)

    ax_list = [carr_ax, low_lat_ax, high_lat_ax, square_wind_ax, dot_wind_ax, hist_ax]

    import string
    # Iterate over the axes objects and their indices
    for ii, this_ax in enumerate(ax_list):
        # Generate label (a), (b), (c), etc.
        label = f"({string.ascii_lowercase[ii]})"

        # Add text to the axes
        # You can adjust the x, y coordinates and other properties as needed

        yy = 0.8 if ii > 0 else 0.6
        this_ax.text(0.02, yy, label, transform=this_ax.transAxes, fontsize=12, va='bottom', ha='left', zorder=1000)



    for jj in [1,2,3]:
        ax_list[jj].set_xticklabels([])

    # Add a colorbar
    cbar_ax = fig.add_axes([0.91, 0.05, 0.03, 0.32])

    # plotobjs = []
    plotobjs = [ cont1]
    # # plotobjs = [scat1, cont1]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend='both', cmap=cmap, extendfrac=0.1, aspect=15)

    cbar.set_label("Wind Speed [km/s]", labelpad=-53)

    fig.set_size_inches((6.5, 12.5))

    plt.subplots_adjust(
        top=0.940,
        bottom=0.065,
        left=0.12,
        right=0.865,
        hspace=0.360,
        wspace=0)


    # Set the output file names
    method = configs.get("flow_method")
    filename = f"png_cr{CR}_f{nwant}_op{n_open}_radial_wind_{method}.png"
    main_file = f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/{filename}'
    wind_file = f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/cr{CR}_f{nwant}_radial_wind_{method}.npy'
    label_file = f'{dat_dir}/batches/{batch}/data/cr{CR}/wind/np_labels.txt'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(wind_file)):
        os.makedirs(os.path.dirname(wind_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))

    # Save the Figures
    print("\n\t\tSaving figures to disk...", end="")
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace("png", "pdf")
    # plt.show()
    # plt.show()
    plt.savefig(outer_file, dpi=200)
    plt.savefig(outer_pdf, format="pdf", dpi=200)

    plt.close(fig)
    print("Success!")

    print("\nDone with wind plotting!\n")
    print("```````````````````````````````\n\n\n")


    ###################################

    vel0_clean = vel0[inds_good]
    ph0_clean = ph0[inds_good]
    th0_clean = th0[inds_good]
    fr0_clean = fr0[inds_good]
    fr1_clean = fr1[inds_good]
    polarity_clean = polarity[inds_good]

    to_save = [ph0_clean, th0_clean, fr0_clean, vel0_clean,
               ph1_clean, th1_clean, fr1_clean, vel1_clean, polarity_clean]

    np.save(wind_file, np.vstack(to_save))

    with open(label_file, 'w') as file:
        file.write("""ph0_clean, th0_clean, fr0_clean, vel0_clean, ph1_clean, th1_clean, fr1_clean, vel1_clean, polarity_clean\n
                   <------------ Solar Surface --------------> <----------- Outer Boundary ---------------->""")
    # np.save(wind_file.replace(".npy", "_interp.npy"), grid_z1)





########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    # Create the argument parser
    configs = configurations()
    do_one = True
    for rotation in configs["rotations"]:
        parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
        parser.add_argument('--cr',     type=int, default=rotation, help='Carrington Rotation')
        parser.add_argument('--interp', type=str, default="linear")
        parser.add_argument('--show',   type=int, default=configs["verbose"])
        parser.add_argument('--file',   type=str, default=None, help='select the file name')
        parser.add_argument('--dat_dir',type=str, default=configs['data_dir'], help='data directory')
        parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='number of fluxons')
        parser.add_argument('--batch',  type=str, default=configs["batch_name"], help='select the batch name')
        parser.add_argument('--adapt',  type=int, default=configs["adapt"], help='Use ADAPT magnetograms')
        args = parser.parse_args()
        configs = configurations(args=args)

        try:
            plot_wind_map_latitude(configs)
        except FileNotFoundError as e:
            print(e, f"{rotation} failed!")
            pass
        if do_one:
            break
