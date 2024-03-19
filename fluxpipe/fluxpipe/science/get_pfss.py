"""
Module for Generating a Fluxon Mapping from Input GONG-Sourced PFSS Coronal Field Solution
==========================================================================================

This module automates the process of generating a fluxon mapping from a given
GONG-sourced PFSS (Potential-Field Source-Surface) coronal field solution. It
provides options to specify the Carrington Rotation (CR), batch name, reduction
factor, data directory, magnetogram file, and other parameters.

Attributes
----------
None

Functions
---------
get_pfss : function
    Main function that loads the fits file, performs PFSS mapping, gets fluxon locations,
    and traces PFSS field lines.

Examples
--------
Run the script with the following command:
    python get_pfss.py --cr 2220 --batch 'my_batch' --reduce 4 --datdir '/path/to/data' --force 1

Authors
-------
Gilly <gilly@swri.org> (and others!)

Dependencies
------------
os, argparse, numpy, pfss_funcs, pipe_helper
"""

# Import required modules
import os
import argparse
import numpy as np
from fluxpipe.science.pfss_funcs import trace_lines, trace_lines_parallel, load_pfss, compute_pfss, load_and_condition_fits_file, get_fluxon_locations
from fluxpipe.helpers.pipe_helper import shorten_path, configurations
import astropy.constants as const
import astropy.units as u
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord

def get_pfss(configs=None):
    """
    Main function to generate a fluxon mapping from a given GONG-sourced PFSS coronal field solution.

    Parameters
    ----------
    configs : dict, optional
        Configuration parameters. If not provided, defaults are loaded from 'configurations'.

    Returns
    -------
    bool
        Returns True upon successful completion.
    """

    configs = configs or configurations()

    from fluxpipe.helpers.pipe_helper import update_magdir_paths
    configs = update_magdir_paths(configs)

    # Extract arguments or use defaults from configs
    cr =        configs.get("cr")
    nwant =     configs.get("nwant")
    magpath =   configs.get("magpath").format(cr)
    force_trace=configs.get("force", 0)
    datdir =    configs.get("data_dir")
    mag_reduce =configs.get("mag_reduce")
    adapt_select =configs.get("adapt_select")
    batch =     configs.get("batch_name")
    adapt =     configs.get("adapt") == 1
    # print(configs['cr'])

    # print("\n\n\n")

    # Print initial message
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("(py) Running PFSS Code to Trace Footpoints into the Corona")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush=True)

    elapsed = 0

    # Load the fits file and format the data and header
    br_safe, fits_path = load_and_condition_fits_file(magpath, datdir, adapt)



    ###############################################################################
    # Do the PFSS mapping
    from fluxpipe.science.pfss_funcs import load_pfss, compute_pfss

    # Get the fluxon locations
    if configs.get("adapt"):
        inst = "adapt"
        adapt = True
        mag_reduce = "f"+str(adapt_select)
    else:
        inst = "hmi"

    pickle_dir = os.path.join(datdir, "pfss")
    if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
    pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{mag_reduce}_{inst}.pkl")
    pfss_output = load_pfss(pickle_path)
    if not pfss_output:
        pfss_output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss

    floc_dir = f"{datdir}/batches/{batch}/data/cr{cr}/floc"
    floc_path   = f"{floc_dir}/floc_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    open_path   = f"{floc_dir}/floc_open_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    closed_path = f"{floc_dir}/floc_closed_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"

    # print(closed_path, "\n\n\n")
    f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch, configs=configs)

    # Trace pfss field lines
    skip_num = 'x'
    timeout_num = 'x'
    # print("Open Path = ", open_path)
    # print("Closed Path = ", closed_path)

    print("\n\tTracing Open and Closed Fluxons...", end="")
    need = not os.path.exists(open_path) or not os.path.exists(closed_path) or force_trace
    print(f"\n\t\tNeed to trace? {need==True}: ", end='')
    if not need:
        try:
            # Just load the lines
            fl_open = np.loadtxt(open_path)
            fl_closed = np.loadtxt(closed_path)
            flnum_open = len(np.unique(fl_open[:, 0]))+1
            flnum_closed = 2*len(np.unique(fl_closed[:, 0]))
            print("Fluxon Loc data files already exist!")
            print(f"\t\t\t{shorten_path(open_path)}")
            print(f"\t\t\t{shorten_path(closed_path)}")
            print(f"\t\tFootpoints:\t Open: {flnum_open}, Closed: \
                    {flnum_closed}, Total: {flnum_open+flnum_closed}")
        except FileNotFoundError:
            need = True

    if need:
        # Trace the lines now
        fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed  = trace_lines(
                                    pfss_output, (f_lon, f_lat, f_sgn), open_path, closed_path, adapt)


    # Record stats in the pfss_output file
    shp = br_safe.data.shape
    pix = shp[0]*shp[1]
    timefile = f'{datdir}/batches/{batch}/pipe_log.txt'
    with open(timefile, 'a+', encoding="utf-8") as f:
        elap = f"\ncr: {cr}, r: {mag_reduce}, rx: {shp[0]}, ry: {shp[1]}, \
                pf_elp: {elapsed:0>3.3f}, t_kpix: {1000*elapsed/pix:0.3f}"
        nlines = f"TrOpen: {flnum_open}, TrClosed: {flnum_closed}, TrGood: \
                {flnum_open+flnum_closed}, TrFail: {skip_num+timeout_num}, "
        f.write(f"{elap}, {nlines}")

    print("\n\t\t\t```````````````````````````````\n \n")
    return True

def get_regular_pfss(configs=None):
    """
    Main function to generate a fluxon mapping from a given GONG-sourced PFSS coronal field solution.

    Parameters
    ----------
    configs : dict, optional
        Configuration parameters. If not provided, defaults are loaded from 'configurations'.

    Returns
    -------
    bool
        Returns True upon successful completion.
    """

    configs = configs or configurations()

    from fluxpipe.helpers.pipe_helper import update_magdir_paths
    configs = update_magdir_paths(configs)

    # Extract arguments or use defaults from configs
    cr =        configs.get("cr")
    nwant =     configs.get("nwant")
    magpath =   configs.get("magpath").format(cr)
    force_trace=configs.get("force", 0)
    datdir =    configs.get("data_dir")
    mag_reduce =configs.get("mag_reduce")
    adapt_select =configs.get("adapt_select")
    batch =     configs.get("batch_name")
    adapt =     configs.get("adapt") == 1
    # print(configs['cr'])

    # print("\n\n\n")

    # Print initial message
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("(py) Running PFSS Code to Trace Footpoints into the Corona")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush=True)

    elapsed = 0

    # Load the fits file and format the data and header
    br_safe, fits_path = load_and_condition_fits_file(magpath, datdir, adapt)

    ###############################################################################
    # Do the PFSS mapping
    from fluxpipe.science.pfss_funcs import load_pfss, compute_pfss

    # Get the fluxon locations
    if configs.get("adapt"):
        inst = "adapt"
        adapt = True
        mag_reduce = "f"+str(adapt_select)
    else:
        inst = "hmi"

    pickle_dir = os.path.join(datdir, "pfss")
    if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
    pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{mag_reduce}_{inst}.pkl")
    pfss_output = load_pfss(pickle_path)
    if not pfss_output:
        pfss_output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss
    import numpy as np
    smaller_br = br_safe.data[::8, ::8]
    f_sgn = np.sign(smaller_br)
    import astropy.units as u
    # Create latitude array
    lat_steps = np.linspace(-1,1, f_sgn.shape[0])
    f_lat = np.tile(lat_steps[:, None], (1, f_sgn.shape[1])).flatten() * u.rad

    # Create longitude array
    lon_steps = np.linspace(0, 2*np.pi, f_sgn.shape[1])
    f_lon = np.tile(lon_steps[None, :], (f_sgn.shape[0], 1)).flatten() * u.rad

    f_sgn = f_sgn.flatten()

    # f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch, configs=configs)

    import astropy.constants as const
    import astropy.units as u
    import matplotlib.colors as mcolor
    import matplotlib.pyplot as plt
    import numpy as np
    import sunpy.map
    from astropy.coordinates import SkyCoord

    # Trace pfss field lines
    skip_num = 'x'
    timeout_num = 'x'
    nwant = "Regular"
    floc_dir = f"{datdir}/batches/{batch}/data/cr{cr}/floc_regular"
    if not os.path.exists(floc_dir): os.makedirs(floc_dir)
    open_path   = f"{floc_dir}/floc_open_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    closed_path = f"{floc_dir}/floc_closed_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    need = not os.path.exists(open_path) or not os.path.exists(closed_path) or force_trace

    if not need:
        try:
            # Just load the lines
            fl_open = np.loadtxt(open_path)
            fl_closed = np.loadtxt(closed_path)
            flnum_open = len(np.unique(fl_open[0]))+1
            flnum_closed = 2*len(np.unique(fl_closed[0]))
            print("\t\tFluxon Loc data files already exist!")
            print(f"\t\t\t{shorten_path(open_path)}")
            print(f"\t\t\t{shorten_path(closed_path)}")
            print(f"\t\tFootpoints:\t Open: {flnum_open}, Closed: \
                    {flnum_closed}, Total: {flnum_open+flnum_closed}")
        except FileNotFoundError:
            need = True

    if need:
        # Trace the lines now
        fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed  = trace_lines_parallel(
                                    pfss_output, (f_lon, f_lat, f_sgn), open_path, closed_path, adapt)

    offid, pol, oth, oph, orad = np.loadtxt(open_path).T
    cffid, pol, cth, cph, crad = np.loadtxt(closed_path).T

    # I want the next line to get the unique values of ffid and then give me the indices of the first element for each unique value
    # This will give me the indices of the first element of each fluxon
    unique_oinds = np.unique(offid, return_index=True)[1]
    unique_oth, unique_oph = oth[unique_oinds], oph[unique_oinds]

    unique_cinds = np.unique(cffid, return_index=True)[1]
    unique_cth, unique_cph = cth[unique_cinds], cph[unique_cinds]


    # photo_op = (1.01 > orad) & (orad > 0.99)
    # photo_cl = (1.01 > crad) & (crad > 0.99)
    import matplotlib.pyplot as plt
    # plt.scatter(f_lon, f_lat, s=10, c=f_sgn, cmap='bwr')

    plt.scatter(f_lon, f_lat, s=10, c='g')
    plt.scatter(unique_oph/180*np.pi, np.sin(np.deg2rad(unique_oth)), s=10, c='r')
    plt.scatter(unique_cph/180*np.pi, np.sin(np.deg2rad(unique_cth)), s=10, c='b')


    # plt.scatter(cph[photo_cl]/180*np.pi, np.sin(np.deg2rad(cth[photo_cl])), s=10, c='b')
    # plt.scatter(oph[photo_op]/180*np.pi, np.sin(np.deg2rad(oth[photo_op])), s=10, c='r')
    # plt.imshow(smaller_br, origin='lower', cmap='gray', zorder=-1)
    plt.show()


    print("\n\t\t\t```````````````````````````````\n \n")
    return True

def get_regular_pfss2(configs=None):
    """
    Main function to generate a fluxon mapping from a given GONG-sourced PFSS coronal field solution.

    Parameters
    ----------
    configs : dict, optional
        Configuration parameters. If not provided, defaults are loaded from 'configurations'.

    Returns
    -------
    bool
        Returns True upon successful completion.
    """

    configs = configs or configurations()

    from fluxpipe.helpers.pipe_helper import update_magdir_paths
    configs = update_magdir_paths(configs)

    # Extract arguments or use defaults from configs
    cr =        configs.get("cr")
    nwant =     configs.get("nwant")
    magpath =   configs.get("magpath").format(cr)
    force_trace=configs.get("force", 0)
    datdir =    configs.get("data_dir")
    mag_reduce =configs.get("mag_reduce")
    adapt_select =configs.get("adapt_select")
    batch =     configs.get("batch_name")
    adapt =     configs.get("adapt") == 1
    # print(configs['cr'])

    # print("\n\n\n")

    # Print initial message
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("(py) Running PFSS Code to Trace Footpoints into the Corona")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush=True)

    elapsed = 0

    # Load the fits file and format the data and header
    br_safe, fits_path = load_and_condition_fits_file(magpath, datdir, adapt)

    ###############################################################################
    # Do the PFSS mapping
    from fluxpipe.science.pfss_funcs import load_pfss, compute_pfss

    # Get the fluxon locations
    if configs.get("adapt"):
        inst = "adapt"
        adapt = True
        mag_reduce = "f"+str(adapt_select)
    else:
        inst = "hmi"

    pickle_dir = os.path.join(datdir, "pfss")
    if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
    pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{mag_reduce}_{inst}.pkl")
    pfss_output = load_pfss(pickle_path)
    if not pfss_output:
        pfss_output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss



    # Trace pfss field lines
    nwant = "Regular"
    floc_dir = f"{datdir}/batches/{batch}/data/cr{cr}/floc_regular"
    if not os.path.exists(floc_dir): os.makedirs(floc_dir)
    open_path   = f"{floc_dir}/floc_open_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    closed_path = f"{floc_dir}/floc_closed_cr{cr}_r{mag_reduce}_f{nwant}_{inst}.dat"
    need = not os.path.exists(open_path) or not os.path.exists(closed_path) or force_trace

    if not need:
        try:
            # Just load the lines
            fl_open = np.loadtxt(open_path)
            fl_closed = np.loadtxt(closed_path)
            flnum_open = len(np.unique(fl_open[0]))+1
            flnum_closed = 2*len(np.unique(fl_closed[0]))
            print("\t\tFluxon Loc data files already exist!")
            print(f"\t\t\t{shorten_path(open_path)}")
            print(f"\t\t\t{shorten_path(closed_path)}")
            print(f"\t\tFootpoints:\t Open: {flnum_open}, Closed: \
                    {flnum_closed}, Total: {flnum_open+flnum_closed}")
        except FileNotFoundError:
            need = True

    from pfsspy import tracing

    r = const.R_sun
    # Number of steps in cos(latitude)
    nsteps = 45
    lon_1d = np.linspace(0, 2 * np.pi, nsteps * 2 + 1)
    lat_1d = np.arcsin(np.linspace(-1, 1, nsteps + 1))
    lon, lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')
    lon, lat = lon*u.rad, lat*u.rad
    seeds = SkyCoord(lon.ravel(), lat.ravel(), r, frame=pfss_output.coordinate_frame)

    print('Tracing field lines...')
    tracer = tracing.FortranTracer(max_steps=3000)
    field_lines = tracer.trace(seeds, pfss_output)
    print('Finished tracing field lines')

    fig = plt.figure()
    # m = pfss_in.map
    ax = fig.add_subplot(2, 1, 1)
    ax.imshow(br_safe.data, origin="lower", cmap='gray', vmin=-500, vmax=500)
    ax.set_title('Input GONG magnetogram')

    ax = fig.add_subplot(2, 1, 2)
    cmap = mcolor.ListedColormap(['tab:red', 'black', 'tab:blue'])
    norm = mcolor.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ncolors=3)
    pols = field_lines.polarities.reshape(2 * nsteps + 1, nsteps + 1).T
    ax.contourf(np.rad2deg(lon_1d), np.sin(lat_1d), pols, norm=norm, cmap=cmap)
    ax.scatter(np.rad2deg(lon), np.sin(lat), s=15, c='g')


    new_olat = field_lines.open_field_lines.solar_feet.lat.to(u.rad)
    new_olon = field_lines.open_field_lines.solar_feet.lon.to(u.rad)


    closed_fields = []
    for fl in field_lines.field_lines:
        if not fl.is_open:
            try:
                closed_fields.append((fl.solar_footpoint.lat.to(u.rad).value, fl.solar_footpoint.lon.to(u.deg).value))
            except IndexError:
                continue


    new_clat, new_clon = zip(*closed_fields)
    # new_clat, new_clon = new_clat, new_clon
    ax.scatter(np.rad2deg(new_olon), np.sin(new_olat), s=10, c='r')
    ax.scatter((new_clon), np.sin((new_clat)), s=10, c='b')

    ax.set_ylabel('sin(latitude)')

    ax.set_title('Open (blue/red) and closed (black) field')
    ax.set_aspect(0.5 * 360 / 2)

    plt.show()


    # if need:
    #     # Trace the lines now
    #     fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed  = trace_lines(
    #                                 pfss_output, (f_lon, f_lat, f_sgn), open_path, closed_path, adapt)

    # offid, pol, oth, oph, orad = np.loadtxt(open_path).T
    # cffid, pol, cth, cph, crad = np.loadtxt(closed_path).T

    # # I want the next line to get the unique values of ffid and then give me the indices of the first element for each unique value
    # # This will give me the indices of the first element of each fluxon
    # unique_oinds = np.unique(offid, return_index=True)[1]
    # unique_oth, unique_oph = oth[unique_oinds], oph[unique_oinds]

    # unique_cinds = np.unique(cffid, return_index=True)[1]
    # unique_cth, unique_cph = cth[unique_cinds], cph[unique_cinds]


    # # photo_op = (1.01 > orad) & (orad > 0.99)
    # # photo_cl = (1.01 > crad) & (crad > 0.99)
    # import matplotlib.pyplot as plt
    # # plt.scatter(f_lon, f_lat, s=10, c=f_sgn, cmap='bwr')

    # plt.scatter(f_lon, f_lat, s=10, c='g')
    # plt.scatter(unique_oph/180*np.pi, np.sin(np.deg2rad(unique_oth)), s=10, c='r')
    # plt.scatter(unique_cph/180*np.pi, np.sin(np.deg2rad(unique_cth)), s=10, c='b')


    # # plt.scatter(cph[photo_cl]/180*np.pi, np.sin(np.deg2rad(cth[photo_cl])), s=10, c='b')
    # # plt.scatter(oph[photo_op]/180*np.pi, np.sin(np.deg2rad(oth[photo_op])), s=10, c='r')
    # # plt.imshow(smaller_br, origin='lower', cmap='gray', zorder=-1)
    # plt.show()


    print("\n\t\t\t```````````````````````````````\n \n")
    return True

########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a fluxon mapping from a GONG-sourced PFSS coronal field solution.')
    # configs_t = configurations()
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--nwant', type=int, default=3000, help='Number of fluxons wanted')
    parser.add_argument('--magpath', type=str, default=None, help='Magnetogram file')
    parser.add_argument('--force', type=int, default=0, help='Force computation of PFSS mapping')
    parser.add_argument('--adapt', type=int, default=0, help='Use ADAPT magnetograms')
    parser.add_argument('--regular', type=int, default=False, help='use regular footpoints')
    args = parser.parse_args()
    configs = configurations(args=args)

    if args.regular:
        get_regular_pfss2(configs)
    else:
    # Run the main function
        get_pfss(configs)
