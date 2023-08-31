"""
Download a Magnetogram for a Particular Carrington Rotation
===========================================================

This script automates the process of downloading a magnetogram for a specific Carrington Rotation (CR).
It provides options to specify the CR, data directory, reduction factor, and whether to force download the files.

Usage:
    python get_magnetograms.py [--cr CARRINGTON_ROTATION] [--datdir DATA_DIRECTORY]
                               [--reduce REDUCTION_FACTOR] [--do_download DOWNLOAD_FLAG]

Arguments:
    --cr:           The Carrington Rotation for which the magnetogram is to be downloaded. Default is 2219.
    --datdir:       The directory where the downloaded data will be stored. Default is '/Users/cgilbert/vscode/fluxons/fluxon-data'.
    --reduce:       The factor by which the magnetogram is reduced. Default is 5.
    --do_download:  Flag to force download the files. Default is 0 (False).

Functions:
    get_magnetogram_file:  Fetches the magnetogram file based on the provided arguments.
    write_magnetogram_params: Writes the parameters of the downloaded magnetogram to a file.

Example:
    python get_magnetograms.py --cr 2220 --datdir '/path/to/data' --reduce 4 --do_download 1

Author:
    Gilly <gilly@swri.org> (and others!)

Dependencies:
    py_pipe_helper
"""

import argparse
from py_pipe_helper import get_magnetogram_file, write_magnetogram_params

datdir = '/Users/cgilbert/vscode/fluxons/fluxon-data'


# create the argument parser
parser = argparse.ArgumentParser(
    description='This script downloads a magnetogram for a particular Carrington Rotation')
parser.add_argument('--cr', type=int, default=2219, help='Carrington Rotation')
parser.add_argument('--datdir', type=str,
                    default=datdir, help='data directory')
parser.add_argument('--reduce', type=int, default=5,
                    help='factor by which the magnetogram is reduced')
parser.add_argument('--do_download', type=int,
                    default=0, help='download the files')
args = parser.parse_args()

# get the magnetogram files
big_path, small_path = get_magnetogram_file(
    cr=args.cr, datdir=args.datdir, force_download=args.do_download, reduce=args.reduce)

# write the magnetogram parameters
write_magnetogram_params(args.datdir, args.cr, small_path, args.reduce)
