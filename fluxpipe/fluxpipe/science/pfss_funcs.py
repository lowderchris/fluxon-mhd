"""
Primary Library for PFSS (Potential Field Source Surface) Modeling.
=====================================================================

This module provides utility functions for:

- Loading and conditioning FITS files containing magnetogram data.
- Converting pixel coordinates to latitude and longitude.
- Computing the PFSS solution based on the conditioned magnetogram.
- Tracing magnetic field lines.
- Plotting fluxon locations on the magnetogram.
- Saving and loading PFSS solutions from pickle files.

Attributes:
    None

Functions:
    load_fits(file_path: str) -> object
    condition_data(data: object) -> object
    compute_pfss(data: object) -> object
    trace_field_lines(pfss_solution: object) -> list
    plot_fluxons(magnetogram: object, field_lines: list) -> None
    save_pfss_solution(pfss_solution: object, file_path: str) -> None
    load_pfss_solution(file_path: str) -> object

Dependencies:
    - astropy
    - matplotlib
    - numpy
    - pfsspy
    - sunpy
    - timeout_decorator
    - tqdm


Author:

    Gilly <gilly@swri.org> (and others!)

"""


import copy
import os
import pickle
from os import path
from time import time

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pfsspy
import sunpy.map
import timeout_decorator
from pfsspy import tracing
from tqdm import tqdm

from fluxpipe.helpers.pipe_helper import load_fits_magnetogram, read_fits_data, shorten_path

mpl.use("qt5agg")

# A note from the original author:
    # We can now use SunPy to load the HMI fits file, and extract the magnetic field data. Interpolate the data down to a more reasonable resolution.
    # The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is
    # not obvious this is the correct way to do this, so use the following lines
    # at your own risk!
    # # brmap = (sunpy.io.fits.read(datdir + fname)) # b radial map
    # # brdat = brmap[0].data[2,:,:] # b data itself



def load_and_condition_fits_file(fname, datdir, adapt):
    """
    Load and condition the FITS file containing the magnetogram data.
    This involves subtracting the mean from the data to enforce div(B) = 0 on the solar surface.
    The metadata is updated to reflect the correct coordinate system, the nans are removed,
    and the data is returned as a sunpy map object.

    Parameters
    ----------
    fname : str
        Filename path
    datdir : str
        Data directory path

    Returns
    -------
    br_safe: GenericMap
        Processed magnetogram data
    fits_path: str
        FITS file path
    """

    # Load the fits file and format the data and header
    print("\n\tLoading and conditioning fits file...")

    try:
        fits_path = fname
        hdulist = read_fits_data(fits_path)
    except FileNotFoundError:
        fits_path = os.path.join(datdir, fname)
        hdulist = read_fits_data(fits_path)
    print("\t\t", shorten_path(fits_path))

    hdu_0 = hdulist[0]
    brdat = hdu_0.data
    header = hdu_0.header


    brdat = brdat - np.mean(brdat)
    radial_magnetogram = sunpy.map.Map(brdat, header)

    radial_magnetogram.meta['ctype1'] = 'CRLN-CEA'
    radial_magnetogram.meta['ctype2'] = 'CRLT-CEA'
    radial_magnetogram.meta['naxis'] = 2

    radial_magnetogram.meta['crlt_obs'] = 0.0
    radial_magnetogram.meta['crln_obs'] = 0.0
    radial_magnetogram.meta['dsun_obs'] = 1.0

    br_safe = copy.copy(radial_magnetogram)
    br_safe.data[np.isnan(br_safe.data)] = 0

    return br_safe, fits_path


def fits_path_to_pickle_path(fits_path, reduce):
    """Returns the path to the pickle file, given a fits file path

    Parameters
    ----------
    fits_path : str
        FITS file path
    reduce : int
        Reduction factor

    Returns
    -------
    str
        Output file path
    """

    the_dir = path.dirname(path.dirname(fits_path))
    pickle_path = path.join(the_dir, f"pfss_output_{reduce}.pkl")
    return pickle_path


def load_pfss(pickle_path):
    """Loads the pickle file

    Parameters
    ----------
    pickle_path : str
        PFSS solution saved state path

    Returns
    -------
    output
        PFSS solution state function
    """

    print("\n\tGetting PFSS...", end="")
    try:
        with open(pickle_path, 'rb') as inp:
            output = pickle.load(inp)
            print(f"Success! Loaded: \n\t\t {shorten_path(pickle_path)}")
        return output
    except FileNotFoundError:
        print("File not found.")
        return None


def compute_pfss(br_safe, pickle_path, nrho=60, rss=2.5):
    """The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
        rho = ln(r), and r is the standard spherical radial coordinate. We need to
        define the number of rho grid points, and the source surface radius.

    Parameters
    ----------
    br_safe : HDU
        Magnetogram hdu
    pickle_path : str
        PFSS solution saved state path
    nrho : int, optional
        Number of points in the rho direction. Defaults to 60.
    rss : float, optional
        Source surface height in r/r_sun. Defaults to 2.5.

    Returns
    -------
    output
        PFSS solution state function
    elapsed : float
        Elapsed runtime in seconds
    """

    # Ingest inputs
    elapsed = 0
    shp = br_safe.data.shape
    print(f"\t\tComputing PFSS on {shp} magnetogram (fairly slow)...", end="", flush=True)
    before = time()

    ###############################################################################
    # From the boundary condition, number of radial grid points, and source
    # surface, we now construct an Input object that stores this information
    pfss_in = pfsspy.Input(br_safe, nrho, rss)

    ###############################################################################
    # Now calculate the PFSS solution
    output = pfsspy.pfss(pfss_in)

    ###############################################################################
    # Report the time taken
    elapsed = time() - before
    print(f"Done! Took {elapsed:.2f} seconds.")

    ###############################################################################
    # Save the results
    with open(pickle_path, 'wb') as outp:
        print("\t\tSaving Pickled Results...", end="")
        pickle.dump(output, outp, pickle.HIGHEST_PROTOCOL)
        print("Success!")
        # print(pickle_path)
    return output, elapsed


def pixel_to_latlon(mag_map, header, fluxon_location):
    """Converts magnetogram pixel index to lat/lon

    Parameters
    ----------
    mag_map : np.ndarray
        Magnetogram data array
    header : Header
        Magnetogram header object
    fluxon_location : Tuple
        Indices of the fluxon root locations

    Returns
    -------
    f_lat : float
        Fluxon root latitudes
    f_lon : float
        Fluxon root longitudes
    f_sgn : float
        Fluxon magnetic sign
    n_flux : int
        Fluxon count
    """

    radial_magnetogram = sunpy.map.Map(mag_map, header)
    shp = mag_map.shape
    n_lat = shp[0]
    lat_center = n_lat // 2
    f_lon = np.deg2rad(fluxon_location[:, 0]
                       * radial_magnetogram.meta['cdelt1'])
    f_lat = (fluxon_location[:, 1]-lat_center) * \
        radial_magnetogram.meta['cdelt2']
    f_sgn = fluxon_location[:, 2]
    n_flux = len(f_sgn)
    return f_lat, f_lon, f_sgn, n_flux


def get_fluxon_locations(floc_path, batch, configs=None):
    """Loads the fluxon location file

    Parameters
    ----------
    floc_path : str
        Path to the fluxon location file
    batch : str
        Name of the current batch

    Returns
    -------
    f_lat : float
        Fluxon root latitudes
    f_lon : float
        Fluxon root longitudes
    f_sgn : float
        Fluxon magnetic sign
    n_flux : int
        Fluxon count
    """

    fluxon_location = np.genfromtxt(floc_path)
    magnet, header = load_fits_magnetogram(batch=batch, ret_all=True, configs=configs)
    f_lat, f_lon, f_sgn, n_flux = pixel_to_latlon(
        magnet, header, fluxon_location)
    return f_lat, f_lon, f_sgn, n_flux


def plot_fluxon_locations(br_safe, cr, datdir, fits_path, reduce,
                          force_plot=False, batch='fluxon', nwant=0, do_plot=False):
    """Plot computed fluxon root locations

    Parameters
    ----------
    br_safe : GenericMap
        Magneogram map object
    cr : str
        Carrington rotation number
    datdir : str
        Data path directory
    fits_path : str
        FITS file path directory
    reduce : int
        Reduction factor
    force_plot : bool, optional
        Plot overwrite toggle. Defaults to False.
    batch : str, optional
        Output file descriptor label. Defaults to 'fluxon'.
    nwant : int, optional
        Output file descriptor label. Defaults to 0.
    do_plot : bool, optional
        Plot toggle. Defaults to False.

    Returns
    -------
    f_lat : float
        Fluxon root latitudes
    f_lon : float
        Fluxon root longitudes
    f_sgn : float
        Fluxon magnetic sign
    n_flux : int
        Fluxon count
    """

    floc_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_cr{cr}_r{reduce}_f{nwant}.dat"
    f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch, cr)
    fluxon_location = np.genfromtxt(floc_path)

    if not do_plot:
        return f_lat, f_lon, f_sgn, n_flux

    the_dir = path.dirname(path.dirname(fits_path))

    top_dir = path.join(datdir, f"batches/{batch}/imgs/footpoints")
    fluxon_map_output_path_top = path.join(
        top_dir, f'cr{cr}_footpoints_r{reduce}_f{nwant}.pdf')

    fluxon_map_output_path_blank = path.join(
        the_dir, "magnetograms", f'footpoints_cr{cr}_r{reduce}_blank.png')
    fluxon_map_output_path_blank_top = path.join(
        top_dir, f'cr{cr}_footpoints_r{reduce}_blank.png')

    if not path.exists(top_dir):
        os.makedirs(top_dir)

    plot_paths = [
        fluxon_map_output_path_top,
        fluxon_map_output_path_blank_top,
        fluxon_map_output_path_blank,
    ]

    print("\n\tPlotting Fluxon Footpoint Locations...")

    need_plot = False
    for testpath in plot_paths:
        if not path.exists(testpath):
            need_plot = True
    if need_plot or force_plot:
        # Print the Fluxon Map
        fig, ax = plt.subplots()

        magnet = br_safe.data
        # find the max and min of the magnetogram plot for use in setting the colormap,
        sigma = 2
        mmean = np.nanmean(magnet)
        msig = np.nanstd(magnet)

        mvmin = mmean - sigma*msig
        mvmax = mmean + sigma*msig
        mvmin, mvmax = -500, 500

        magimg1 = ax.imshow(magnet, cmap='gray', interpolation=None,
                            origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)
        shp = magnet.shape

        plt.axis('off')
        sz0 = 6  # inches
        ratio = shp[1]/shp[0]
        sz1 = sz0*ratio  # inches
        DPI = shp[1] / sz1  # pixels/inch
        fig.set_size_inches((sz1, sz0))
        plt.tight_layout()

        # plot a blank version of the map
        plt.savefig(fluxon_map_output_path_blank,
                    bbox_inches='tight', dpi=4*DPI)

        # plot a blank version of the map
        magimg2 = ax.imshow(magnet, cmap='gray', interpolation=None,
                            origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)
        plt.savefig(fluxon_map_output_path_blank_top,
                    bbox_inches='tight', dpi=4*DPI)

        # scatter the fluxons on top of the map
        magimg3 = ax.imshow(magnet, cmap='gray', interpolation=None,
                            origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)

        x, y, sig = zip(*fluxon_location)
        colors = ['red' if s > 0 else 'teal' for s in sig]
        ax.scatter(x, y, c=colors, alpha=0.4)

        # plot the fluxons scattered on top of the map
        plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)
        plt.close()
        print("Success! Files saved to: ", end="\n\t\t")
    else:
        print("Skipped! Files already exist:", end="\n\t\t")

    for testpath in plot_paths:
        print(shorten_path(testpath), end="\n\t\t")

    return f_lat, f_lon, f_sgn, n_flux


def trace_lines(output, f_vars, open_path, closed_path, adapt):
    """ This function traces field lines using the PFSS solution.

    Parameters
    ----------
    output
        PFSS solution state function
    f_vars : Tuple
        Storage of fluxon longitude, latitude, and magnetic sign
    open_path : str
        File path for writing open fluxon coordinates
    closed_path : str
        File path for writing closed fluxon coordiantes

    Returns
    -------
    fl_open : np.ndarray
        Open feldline coordinates
    fl_closed : np.ndarray
        Closed fieldline coordiantes
    skip_num : int
        Counter of skipped tracings
    timeout_num : int
        Counter of timedout tracings
    flnum_open : int
        Count of open fieldlines
    flnum_closed : int
        Count of closed fieldlines
    """

    # TODO Note that this code was originally developed using an older version of the pfsspy code -
    # improvements look to be able to handle bundles of fieldlines,
    # which would greatly simplify some of the code. Do this.
    (f_lon, f_lat, f_sgn) = f_vars

    # Capture an array of starting and endpoints for open fieldlines 2x(r, lat, lon)
    # + sgn of starting point
    fl_open = np.zeros([1, 5])
    flnum_open = 0

    # Capture an array of starting and endpoints for closed fieldlines 2x(r, lat, lon)
    # + sgn of starting point
    fl_closed = np.zeros([1, 5])
    flnum_closed = 0

    r0 = 1.01 * const.R_sun

    skip_num = 0
    timeout_num = 0
    tracer = tracing.PythonTracer()
    for i, coords in enumerate(tqdm(zip(f_lon, f_lat),
                            desc="Tracing Field Lines", total=len(f_lat))):
        try:
            output, fl_open, fl_closed, flnum_open, flnum_closed, tracer = trace_each(
                coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer, r0, f_sgn, adapt)
        except timeout_decorator.TimeoutError:
            timeout_num += 1
        except ValueError:
            skip_num += 1

    if skip_num > 0 or timeout_num > 0:
        t_perc = 100*timeout_num/len(f_lon)
        s_perc = 100*skip_num/len(f_lon)
        print(f"\n\t\tSome iterations failed. Timed-out: {timeout_num} \
                ({t_perc:0.2f}%), ValueError: {skip_num} ({s_perc:0.2f}%)\n")

    print(f"\nOpen Lines: {flnum_open+1}, Closed Lines: {flnum_closed}, \
          Failures: {skip_num+timeout_num}, Total Good: {flnum_open+flnum_closed}")

    fl_open = fl_open[1:]
    fl_closed = fl_closed[1:]

    # Save these coordinates to file
    print("\n\tSaving Fluxons...", end="")
    np.savetxt(open_path, fl_open)
    np.savetxt(closed_path, fl_closed)
    print("Success!")

    return fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed



@timeout_decorator.timeout(5)
def trace_each(coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer, r0, f_sgn, adapt):
    # import pdb; pdb.set_trace()
    """Fieldline tracer for each fieldline

    Parameters
    ----------
    coords : Tuple
        Set of input fieldline longitude and latitudes
    i : int
        Fieldline tracing counter
    output
        PFSS solution state function
    fl_open : np.ndarray
        Open feldline coordinates
    fl_closed : np.ndarray
        Closed fieldline coordiantes
    flnum_open : int
        Count of open fieldlines
    flnum_closed : int
        Count of closed fieldlines
    tracer
        PFSSPy tracing function
    r0 : float
        Starting radius
    f_sgn : float
        Fluxon magnetic sign

    Returns
    -------
    output
        PFSS solution state function
    fl_open : np.ndarray
        Open feldline coordinates
    fl_closed : np.ndarray
        Closed fieldline coordiantes
    flnum_open : int
        Count of open fieldlines
    flnum_closed : int
        Count of closed fieldlines
    tracer
        PFSSPy tracing function
    """

    (this_flon, this_flat) = coords
    if adapt:
        # this_flon = np.deg2rad(this_flon)
        this_flat = np.deg2rad(this_flat)
    # import pdb; pdb.set_trace()
    coord_frame = output.coordinate_frame
    # THE ARCSIN IS A TEST
    x0 = SkyCoord(this_flon * u.rad, np.arcsin(this_flat)
                  * u.rad, r0, frame=coord_frame)

    fl = output.trace(tracer, x0)
    fl = fl.field_lines[0]
    if fl.is_open:
        fl_pol = fl.polarity
        fl = fl.coords.spherical
        fl_lats = fl.lat.value
        fl_lons = fl.lon.value
        fl_rads = (fl.distance/const.R_sun).value

        if fl_pol == -1:
            (fl_lats, fl_lons, fl_rads) = (
                fl_lats[::-1], fl_lons[::-1], fl_rads[::-1])
        fl_lats = np.concatenate((fl_lats, fl_lats[-1]*np.ones(10)), axis=0)
        fl_lons = np.concatenate((fl_lons, fl_lons[-1]*np.ones(10)), axis=0)
        fl_rads = np.concatenate((fl_rads, np.linspace(2.5, 22, 10)), axis=0)
        prev_rad = 0.
        for j in np.arange(0, len(fl_lats)):
            # Output flnum, polarity, latitude, longitude, radius
            if np.abs(fl_rads[j] - prev_rad) > 0.1:
                fl_open = np.append(
                    fl_open, [[flnum_open, f_sgn[i], fl_lats[j], fl_lons[j], fl_rads[j]]], axis=0)
                prev_rad = fl_rads[j]
        flnum_open += 1
    else:
        fl = fl.coords.spherical
        fl_rads = (fl.distance/const.R_sun).value
        max_rad = fl_rads.max()
        rad_thresh = np.max([((max_rad-1.) / 5.), 0.01])
        prev_rad = 0.
        if np.mod(flnum_closed, 2) == 0:
            for j in np.arange(0, len(fl.lat)):
                # Output flnum, polarity, latitude, longitude, radius
                if (np.abs(fl_rads[j] - prev_rad) > rad_thresh) or (fl_rads[j] == max_rad):
                    fl_closed = np.append(fl_closed, [[flnum_closed, f_sgn[i], fl.lat[j].value,
                                    fl.lon[j].value, (fl.distance[j]/const.R_sun).value]], axis=0)
                    prev_rad = fl_rads[j]
        flnum_closed += 1
    return output, fl_open, fl_closed, flnum_open, flnum_closed, tracer
