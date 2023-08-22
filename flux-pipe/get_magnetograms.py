# Module Docstring
"""Download a magnetogram for a particular Carrington Rotation.
Args:
    cr (int): Carrington Rotation
    datdir (str): data directory
    reduce (int): factor by which the magnetogram is reduced
    do_download (int): download the files
Returns:
    None
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
