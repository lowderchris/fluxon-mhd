"""
Changed
==========================================================

This script is designed to plot the magnetogram along with footpoints. It provides
options to specify the Carrington Rotation (CR), batch name, reduction factor, data directory,
and other parameters.

Usage:
    python plot_fieldmap.py [--cr CARRINGTON_ROTATION] [--nwant NUMBER_WANTED]
                             [--open OPEN_FLUXONS] [--closed CLOSED_FLUXONS]
                             [--dat_dir DATA_DIRECTORY] [--batch BATCH_NAME]

Arguments:
    --cr:           The Carrington Rotation for which the magnetogram is to be plotted. Default is None.
    --nwant:        The number of fluxons wanted. Default is None.
    --open:         The number of open fluxons. Default is None.
    --closed:       The number of closed fluxons. Default is None.
    --dat_dir:      The directory where the data will be stored. Default is defined in the config.ini file.
    --batch:        The batch name for the operation. Default is 'default_batch'.

Functions:
    magnet_plot:    A primary function for plotting the magnetogram with footpoints.

Example:
    python plot_fieldmap.py --cr 2220 --nwant 100 --open 50 --closed 50 --dat_dir '/path/to/data' --batch 'my_batch'

Author:
    Gilly <gilly@swri.org> (and others!)

Dependencies:
    os, os.path, argparse, matplotlib.pyplot, numpy, pfss_funcs, pipe_helper
"""

import os
import os.path as path
import argparse
# import matplotlib as mpl; mpl.use("qt5agg")
import matplotlib.pyplot as plt
import numpy as np


from fluxpipe.science.pfss_funcs import pixel_to_latlon
from fluxpipe.helpers.pipe_helper import (configurations, load_fits_magnetogram, load_magnetogram_params,
                            shorten_path, get_ax)

def magnet_plot(
        get_cr=None, datdir=None, _batch=None, open_f=None, closed_f=None, force=False, reduce_amt=0,
        nact=0, nwant=None, do_print_top=True, ax=None, verb=True, ext="png",
        plot_all=True, plot_open=True, do_print=False, vmin=-500, vmax=500, configs=None, legend=False):

    """ This function has been re-imagined into a

    Parameters
    ----------
    get_cr : int
        the carrington rotation number


    Returns
    -------

    """

    fig, ax0 = get_ax(ax)

    if True:
        print("\t\t(py) Determining Footpoint Distances")



    if configs is not None:
        datdir = datdir or configs.get("data_dir", None)
        _batch = _batch or configs.get("batch_name", None)
        get_cr = get_cr or configs.get("cr", None)
        nwant = nwant   or configs.get("nwant", None)
        reduce_amt = reduce_amt or configs.get("mag_reduce", None)
        if configs['adapt']:
            inst = "adapt"
            reduce_amt = "f" + str(configs.get("adapt_select"))
        else:
            inst = "hmi"
    else:
        print("No configs given!")
        raise ValueError


    # Define the directory paths for the files
    floc_path = f"{datdir}/batches/{_batch}/data/cr{get_cr}/floc/"
    top_dir   = f"{datdir}/batches/{_batch}/imgs/footpoints/"
    if not path.exists(top_dir):
        os.makedirs(top_dir)


    # Define the file names with their complete paths
    open_file   = open_f     or   f"{floc_path}floc_open_cr{get_cr}_r{reduce_amt}_f{nwant}_{inst}.dat"
    closed_file = closed_f   or   f"{floc_path}floc_closed_cr{get_cr}_r{reduce_amt}_f{nwant}_{inst}.dat"
    magnet_file = f"{datdir}/magnetograms/CR{get_cr}_r{reduce_amt}_{inst}.fits"
    all_file    = closed_file.replace("closed_", "")
    fname = magnet_file

    # Load the data
    if do_print_top:
        print(f"\t\tOpening {shorten_path(all_file)}...")
    fluxon_location = np.genfromtxt(all_file)
    magnet, header = load_fits_magnetogram(batch=_batch, ret_all=True, configs=configs, fname=fname)
    f_lat, f_lon, f_sgn, _fnum = pixel_to_latlon(magnet, header, fluxon_location)


    if do_print_top:
        print(f"\t\tOpening {shorten_path(open_file)}...")
    oflnum, oflx, olat, olon, orad = np.loadtxt(open_file, unpack=True)

    if do_print_top:
        print(f"\t\tOpening {shorten_path(closed_file)}...\n")
    cflnum, cflx, clat, clon, crad = np.loadtxt(closed_file, unpack=True)


    ## Keep only the values where the radius is 1.0
    rtol = 0.001
    get_r = 1.0

    #Open fields
    oflnum_low = oflnum[np.isclose(orad, get_r, rtol)]
    oflx_low =     oflx[np.isclose(orad, get_r, rtol)]
    olat_low =     olat[np.isclose(orad, get_r, rtol)]
    olon_low =     olon[np.isclose(orad, get_r, rtol)]

    # Closed fields
    cflnum_low = cflnum[np.isclose(crad, get_r, rtol)]
    cflx_low =     cflx[np.isclose(crad, get_r, rtol)]
    clat_low =     clat[np.isclose(crad, get_r, rtol)]
    clon_low =     clon[np.isclose(crad, get_r, rtol)]

    # Convert to radians
    ph_olow, th_olow = np.sin(np.deg2rad(olat_low)), np.deg2rad(olon_low)
    ph_clow, th_clow = np.sin(np.deg2rad(clat_low)), np.deg2rad(clon_low)

    # Report the number of open and closed fluxons
    _n_open = int(np.max(oflnum_low))
    _n_closed = int(np.max(cflnum_low))
    _n_flux = _n_open + _n_closed
    _n_outliers = np.abs(_fnum-_n_flux)
    print(f"\t\t\tOpen: {_n_open}, Closed: {_n_closed}, Total: {_n_flux}, Outliers: {_n_outliers}")

    # Define the file name for the plot
    pic_name = f'cr{get_cr}_f{nwant}_ou{_n_open}_footpoints_topology.{ext}'
    fluxon_map_output_path =   path.join(floc_path, pic_name)
    fluxon_map_output_path_top = path.join(top_dir, pic_name)
    fluxon_csv_output_path_top = path.join(floc_path, "distances.csv")


    # Check if the plot already exists
    do_plot = False
    pic_paths = [fluxon_map_output_path, fluxon_map_output_path_top]
    # pic_paths = [fluxon_map_output_path_top]
    for testpath in pic_paths:
        if not path.exists(testpath):
            do_plot = True
            break

    force=False

    if do_print:
        print("\tPlotting...", end="")
    if do_plot or force or (ax is not None):
        # Plot the magnetogram

        # ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
        #         extent=(0,2*np.pi,-1,1), aspect='auto', vmin=vmin, vmax=vmax, zorder=5, alpha=0.8)
        ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
                extent=(0,2*np.pi,-1,1), aspect='auto', vmin=vmin, vmax=vmax, zorder=-5, alpha=1)



        # # Plot all the fluxons
        # Filter positive and negative cases for closed fluxons
        positive_indices_closed = [i for i, s in enumerate(cflx_low) if s > 0]
        negative_indices_closed = [i for i, s in enumerate(cflx_low) if s < 0]

        # Data for positive and negative cases
        f_lon_positive_closed = [th_clow[i] for i in positive_indices_closed]
        f_lat_positive_closed = [ph_clow[i] for i in positive_indices_closed]
        f_lon_negative_closed = [th_clow[i] for i in negative_indices_closed]
        f_lat_negative_closed = [ph_clow[i] for i in negative_indices_closed]

        if plot_all:
            # Plot positive cases with labels
            ax0.scatter(f_lon_positive_closed, f_lat_positive_closed, s=6**2, c='orange', alpha=0.8, label='Positive')

            # Plot negative cases with labels
            ax0.scatter(f_lon_negative_closed, f_lat_negative_closed, s=6**2, c='teal', alpha=0.8, label='Negative')

        # Filter positive and negative cases for open fluxons
        positive_indices_open = [i for i, s in enumerate(oflx_low) if s > 0]
        negative_indices_open = [i for i, s in enumerate(oflx_low) if s <= 0]

        # Data for positive and negative cases
        f_lon_positive_open = [th_olow[i] for i in positive_indices_open]
        f_lat_positive_open = [ph_olow[i] for i in positive_indices_open]
        f_lon_negative_open = [th_olow[i] for i in negative_indices_open]
        f_lat_negative_open = [ph_olow[i] for i in negative_indices_open]

        if plot_open:
            ax0.scatter(f_lon_positive_open, f_lat_positive_open, s=5**2, c='red', alpha=1.0, label='Positive (Open)', edgecolors='k')
            ax0.scatter(f_lon_negative_open, f_lat_negative_open, s=5**2, c='blue', alpha=1.0, label='Negative (Open)', edgecolors='k')


        # fig, ax2 = plt.subplots()
        # make a voronoi map of the closed footpoints
        # from scipy.spatial import Voronoi, voronoi_plot_2d
        # vor = Voronoi(np.array([th_clow, ph_clow]).T)
        # voronoi_plot_2d(vor, ax=ax0, show_vertices=False, line_colors='orange',
        #         line_width=2, line_alpha=0.6, point_size=2, point_alpha=0.6)

        # for each of the voronoi regions, calculate the distance to the point that owns that region
        # and plot it as a color


        # import numpy as np
        # import matplotlib.pyplot as plt
        from scipy.spatial import cKDTree

        # Provided points (longitude, latitude)
        points = np.array([th_clow, ph_clow]).T #These are in radians and sin(radians)

        # Image dimensions
        img_width, img_height = magnet.T.shape

        # Normalize and scale points to image dimensions
        # Longitude: 0 to 2pi maps to 0 to img_width
        # Latitude: -1 to 1 maps to 0 to img_height
        scaled_points = np.empty_like(points)
        scaled_points[:, 0] = (points[:, 0] / (2 * np.pi)) * img_width
        scaled_points[:, 1] = ((points[:, 1] + 1) / 2) * img_height
        #These are now in pixel units.

        # Create an array with the same shape as the image
        distance_array = np.zeros((img_height, img_width), dtype=float)

        # Generate grid points based on the image shape
        y_indices, x_indices = np.indices((img_height, img_width))
        grid_points = np.column_stack((x_indices.ravel(), y_indices.ravel()))

        # Build a KDTree for efficient nearest-neighbor query
        tree = cKDTree(scaled_points)

        # Calculate degrees per pixel
        deg_per_pixel_x = 360 / img_width  # For longitude
        deg_per_pixel_y = 180 / img_height  # For latitude

        # Query the nearest distance for each grid point, requesting separate x and y components
        distances, _ = tree.query(grid_points, k=1, p=2, workers=-1, eps=0)

        # Convert pixel distances to degrees
        distances_x = (distances // img_width) * deg_per_pixel_x
        distances_y = (distances % img_width) * deg_per_pixel_y
        distances_in_degrees = np.sqrt(distances_x**2 + distances_y**2)

        # Reshape and display
        distance_array_degrees = distances_in_degrees.reshape((img_height, img_width))
        im = ax0.imshow(distance_array_degrees, cmap='viridis', origin='lower', alpha=0.75, extent=(0, 2*np.pi, -1, 1), aspect='auto')
        plt.colorbar(im, label='Distance to Nearest Closed Footpoint (Degrees)')
        plt.title('Distance to Nearest Closed Footpoint (Degrees)')


        from scipy import interpolate

        # Assuming distance_array_degrees is already computed as in the previous steps
        # and img_width, img_height are defined as before

        # Create a grid of longitude and latitude values corresponding to the image
        latitude = np.linspace(-1, 1, img_height)  # -90 to 90 degrees
        longitude = np.linspace(0, 2*np.pi, img_width)  # 0 to 360 degrees

        # Create the interpolator function
        dist_interp = interpolate.RectBivariateSpline(latitude, longitude, distance_array_degrees)
        open_points = np.array([th_olow, ph_olow]).T
        #use the interpolator to get the value at each point
        distances = dist_interp.ev(open_points[:, 1], open_points[:, 0])
        ax0.scatter(th_olow, ph_olow, c=distances, s=10**2, cmap='viridis', alpha=1, edgecolors='r', zorder=100000)

        point_distances = [(th, ph, dist) for th, ph, dist in zip(th_olow, ph_olow, distances)]

        np.savetxt(fluxon_csv_output_path_top, distance_array_degrees, delimiter=", ")

        do_legend=False

        if do_legend:
            ax0.legend(fontsize="small", loc="upper left", framealpha=0.75)

        if ax is None:
            shp = magnet.shape #pixels
            plt.axis('off')
            sz0=6 #inches
            ratio = shp[1]/shp[0]
            sz1=sz0*ratio #inches
            DPI = shp[1] / sz1 #pixels/inch
            fig.set_size_inches((sz1, sz0))
            plt.tight_layout()
            plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)
            # plt.show()
            plt.close(fig)
    else:
        if do_print:
            print("\tSkipped! Files already exist:")
            print(f"\t\t{shorten_path(fluxon_map_output_path)}")
            print(f"\t\t{shorten_path(fluxon_map_output_path_top)}")
    if do_print:
        print(f"\n\t    n_open: {_n_open}, n_closed: {_n_closed}, \
                n_total: {_n_flux}, n_all: {_fnum}, n_outliers: {_n_outliers}")

    if do_print_top:
        print("\t\t    Success!")
        print("\t\t\t```````````````````````````````\n\n")
    return _n_open, _n_closed, _n_flux, _fnum, _n_outliers





########################################################################
# Main Code
# ----------------------------------------------------------------------
#

if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description=
            'This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--file', type=str, default=None, help='Data File Name')
    parser.add_argument('--nwant', type=int, default=None, help='Number of Fluxons')
    parser.add_argument('--open', type=str, default=None)
    parser.add_argument('--closed', type=str, default=None)
    parser.add_argument('--adapt', type=str, default=None)

    args = parser.parse_args()
    configs = configurations(debug=False, args=args)


    # Run the code
    # (hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)

    # magnet_plot(CR,              configs['data_dir'],  configs['batch_name'],
    #             configs['open'], configs['closed'],    do_print=configs['verbose'],
    #             reduce_amt=configs['mag_reduce'],      nwant=nwant,
    #             do_print_top=True)

    magnet_plot(configs=configs)