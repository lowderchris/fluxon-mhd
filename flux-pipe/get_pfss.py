"""
Generate a Fluxon Mapping from Input GONG-Sourced PFSS Coronal Field Solution
=============================================================================

This script automates the process of generating a fluxon mapping from a given
GONG-sourced PFSS (Potential-Field Source-Surface) coronal field solution. It
provides options to specify the Carrington Rotation (CR), batch name, reduction
factor, data directory, magnetogram file, and other parameters.

Usage:
    python get_pfss.py [--cr CARRINGTON_ROTATION] [--batch BATCH_NAME]
                               [--reduce REDUCTION_FACTOR] [--datdir DATA_DIRECTORY]
                               [--magfile MAGNETOGRAM_FILE] [--nwant NUMBER_WANTED]
                               [--force FORCE_FLAG]

Arguments:
    --cr:           The Carrington Rotation for which the PFSS mapping is to be generated. Default is 2219.
    --batch:        The batch name for the operation. Default is 'default_batch'.
    --reduce:       The factor by which the magnetogram is reduced. Default is 5.
    --datdir:       The directory where the data will be stored. Default is None.
    --magfile:      The magnetogram file to be used. Default is None.
    --nwant:        The number of fluxons wanted. Default is None.
    --force:        Flag to force the computation of the PFSS mapping. Default is 0 (False).

Functions:
    trace_lines, load_pfss, compute_pfss, load_and_condition_fits_file, get_fluxon_locations:
                Various utility functions imported from 'pfss_funcs' and 'py_pipe_helper'.

Example:
    python get_pfss.py --cr 2220 --batch 'my_batch' --reduce 4 --datdir '/path/to/data' --force 1

Author:
    Gilly <gilly@swri.org> (and others!)

Dependencies:
    os, argparse, numpy, astropy, pfss_funcs, py_pipe_helper
"""
# First, import required modules
import os
import argparse
import numpy as np
import astropy.constants as const
from pfss_funcs import (trace_lines, load_pfss, compute_pfss, \
                        load_and_condition_fits_file, get_fluxon_locations)
from py_pipe_helper import shorten_path

###############################################################################
# create the argument parser
parser = argparse.ArgumentParser(
    description='This script downloads a magnetogram for a particular Carrington Rotation')
parser.add_argument('--cr', type=int, default=2219, help='Carrington Rotation')
parser.add_argument('--batch', type=str,
                    default='default_batch', help='batch name')
parser.add_argument('--reduce', type=int, default=5,
                    help='factor by which the magnetogram is reduced')
parser.add_argument('--datdir', type=str, default=None, help='data directory')
parser.add_argument('--magfile', type=str, default=None,
                    help='magnetogram file')
parser.add_argument('--nwant', type=int, default=None, help='magnetogram file')
parser.add_argument('--force', type=int, default=0,
                    help='force the computation of the PFSS mapping')

# parse the arguments
args = parser.parse_args()
magfile = args.magfile
datdir = args.datdir
cr = args.cr
reduce = args.reduce
force_plot = force_trace = args.force
batch = args.batch
nwant = args.nwant


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("(py) Running PFSS Code to Trace Footpoints into the Corona")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush=True)

elapsed = 0

###############################################################################
# Load the fits file and format the data and header
br_safe, fits_path = load_and_condition_fits_file(magfile, datdir)


###############################################################################
# Do the PFSS mapping
pickle_dir = os.path.join(datdir, "pfss")
if not os.path.exists(pickle_dir):
    os.makedirs(pickle_dir)
pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{reduce}.pkl")

output = load_pfss(pickle_path)
if not output:
    output, elapsed = compute_pfss(br_safe, pickle_path)  # , nrho, rss)


###############################################################################
# Get the fluxon locations
floc_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_cr{cr}_r{reduce}_f{nwant}.dat"
f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch)


###############################################################################
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


###############################################################################
# Record stats in the output file
shp = br_safe.data.shape
pix = shp[0]*shp[1]
timefile = f'{datdir}/batches/{batch}/pipe_log.txt'
with open(timefile, 'a+', encoding="utf-8") as f:
    # a good name for the variable
    elap = f"\ncr: {cr}, r: {reduce}, rx: {shp[0]}, ry: {shp[1]}, \
            pf_elp: {elapsed:0>3.3f}, t_kpix: {1000*elapsed/pix:0.3f}"
    nlines = f"TrOpen: {flnum_open}, TrClosed: {flnum_closed}, TrGood: \
            {flnum_open+flnum_closed}, TrFail: {skip_num+timeout_num}, "
    f.write(f"{elap}, {nlines}")


print("\n\t\t\t```````````````````````````````\n \n")
