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
from pfss_funcs import (trace_lines, load_pfss, compute_pfss, \
                        load_and_condition_fits_file, get_fluxon_locations)
from pipe_helper import shorten_path, configurations

def get_pfss(configs=None):
    """
    Main function to generate a fluxon mapping from a given GONG-sourced PFSS coronal field solution.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments parsed via argparse.

    Returns
    -------
    bool
        Returns True upon successful completion.
    """
    """Main function to generate a fluxon mapping from a given GONG-sourced PFSS coronal field solution."""
    configs = configs or configurations()

    # Extract arguments or use defaults from configs
    cr =        configs.get("rotations")[0]
    magfile =   configs.get("magfile").format(cr)
    force_trace=configs.get("force", 0)
    datdir =    configs.get("data_dir")
    nwant =     configs.get("fluxon_count")[0]
    reduce =    configs.get("mag_reduce")
    batch =     configs.get("batch_name")

    # Print initial message
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("(py) Running PFSS Code to Trace Footpoints into the Corona")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush=True)

    elapsed = 0

    # Load the fits file and format the data and header
    br_safe, fits_path = load_and_condition_fits_file(magfile, datdir)

    # Do the PFSS mapping
    pickle_dir = os.path.join(datdir, "pfss")
    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir)
    pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{reduce}.pkl")

    output = load_pfss(pickle_path)
    if not output:
        output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss

    # Get the fluxon locations
    floc_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_cr{cr}_r{reduce}_f{nwant}.dat"
    f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch)

    # Trace pfss field lines
    skip_num = 'x'
    timeout_num = 'x'
    open_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_open_cr{cr}_r{reduce}_f{nwant}.dat"
    closed_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_closed_cr{cr}_r{reduce}_f{nwant}.dat"

    print("\n\tTracing Open and Closed Fluxons...", end="")
    if not os.path.exists(open_path) or force_trace:
        trace_out = trace_lines(output, (f_lon, f_lat, f_sgn),
                                open_path, closed_path)
        fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed = trace_out
    else:
        fl_open = np.loadtxt(open_path)
        fl_closed = np.loadtxt(closed_path)
        flnum_open = len(np.unique(fl_open[:, 0]))+1
        flnum_closed = 2*len(np.unique(fl_closed[:, 0]))
        print("Skipped! Floc dat files already exist:")
        print(f"\t\t{shorten_path(open_path)}")
        print(f"\t\t{shorten_path(closed_path)}")
        print(
            f"\t\tFootpoints:\t Open: {flnum_open}, Closed: \
                {flnum_closed}, Total: {flnum_open+flnum_closed}")

    # Record stats in the output file
    shp = br_safe.data.shape
    pix = shp[0]*shp[1]
    timefile = f'{datdir}/batches/{batch}/pipe_log.txt'
    with open(timefile, 'a+', encoding="utf-8") as f:
        elap = f"\ncr: {cr}, r: {reduce}, rx: {shp[0]}, ry: {shp[1]}, \
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
    parser = argparse.ArgumentParser(description='This script downloads a magnetogram for a particular Carrington Rotation')

    # Use default values from configs
    configs = configurations()
    parser.add_argument('--cr',     type=int, default=configs.get('rotations')[0], help='Carrington Rotation')
    parser.add_argument('--nwant',  type=int, default=configs.get('fluxon_count')[0], help='magnetogram file')
    parser.add_argument('--magfile',type=str, default=None, help='magnetogram file')
    parser.add_argument('--force',  type=int, default=0, help='force the computation of the PFSS mapping')
    # parser.add_argument('--batch',  type=str, default=configs.get('batch_name', 'default_batch'), help='batch name')
    # parser.add_argument('--reduce', type=int, default=configs.get('mag_reduce', 2), help='factor by which the magnetogram is reduced')
    # parser.add_argument('--datdir', type=str, default=configs.get('data_dir', None), help='data directory')

    args = parser.parse_args()

    # Update configs with command-line arguments
    for arg, value in vars(args).items():
        if value is not None:
            configs[arg] = value

    # Run the main function
    get_pfss(configs)
