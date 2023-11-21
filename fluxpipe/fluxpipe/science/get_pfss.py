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
from fluxpipe.science.pfss_funcs import trace_lines, load_pfss, compute_pfss, load_and_condition_fits_file, get_fluxon_locations
from fluxpipe.helpers.pipe_helper import shorten_path, configurations

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

    output = load_pfss(pickle_path)
    if not output:
        output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss

    floc_dir = f"{datdir}/batches/{batch}/cr{cr}/floc"
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
                                    output, (f_lon, f_lat, f_sgn), open_path, closed_path, adapt)


    # Record stats in the output file
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




########################################################################
# Main Code
# ----------------------------------------------------------------------
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a fluxon mapping from a GONG-sourced PFSS coronal field solution.')
    # configs_t = configurations()
    parser.add_argument('--cr', type=int, default=2270, help='Carrington Rotation')
    parser.add_argument('--nwant', type=int, default=100, help='Number of fluxons wanted')
    parser.add_argument('--magpath', type=str, default=None, help='Magnetogram file')
    parser.add_argument('--force', type=int, default=0, help='Force computation of PFSS mapping')
    parser.add_argument('--adapt', type=int, default=0, help='Use ADAPT magnetograms')
    args = parser.parse_args()
    configs = configurations(args=args)

    # print("\n\n\n\n\n\n|\n|\n|\n A:")
    # print(args.magpath)
    # print(configs["magpath"])
    # print("\n\n\n\n\n\n|\n|\n|\n")

    # Run the main function
    get_pfss(configs)
