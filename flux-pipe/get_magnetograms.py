"""Download a Magnetogram for a Particular Carrington Rotation
===========================================================
import argparse


This module automates the process of downloading a magnetogram for a specific Carrington Rotation (CR).
It provides options to specify the CR, data directory, reduction factor, and whether to force download the files.

Attributes
----------
configs : dict
    The configuration parameters loaded from 'config.ini'.

Functions
---------
get_magnetogram_file : function
    Fetches the magnetogram file based on the provided arguments.

Examples
--------
Run the script with the following command:
    python get_magnetograms.py --cr 2220  --reduce 4

Authors
-------
Gilly <gilly@swri.org> (and others!)
"""


import argparse
from pipe_helper import (configurations, get_magnetogram_file, get_adapt_file)
configs = configurations()
# Create the argument parser
parser = argparse.ArgumentParser(
    description='This script downloads a magnetogram for a particular Carrington Rotation')
parser.add_argument('--cr', type=int,     default=configs["rotations"][0],  help='Carrington Rotation')
parser.add_argument('--reduce', type=int, default=configs["mag_reduce"],    help='Reduction factor')
args = parser.parse_args()



# get the magnetogram files
get_adapt = True
if get_adapt:
    big_path, small_path = get_ADAPT_file(cr=args.cr, datdir=args.datdir, force_download=args.do_download, reduce=args.reduce)
else:
    big_path, small_path = get_magnetogram_file(cr=args.cr, datdir=args.datdir, force_download=args.do_download, reduce=args.reduce)


