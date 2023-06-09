import argparse
import os
from py_pipe_helper import get_magnetogram_file, reduce_fits_image, find_file_with_string, write_params_file, reduce_mag_file, make_mag_dir
magneto_file = None

# create the argument parser
parser = argparse.ArgumentParser(description='This script downloads a magnetogram for a particular Carrington Rotation')
parser.add_argument('--cr', type=int, default=2219, help='Carrington Rotation')
parser.add_argument('--datdir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
parser.add_argument('--reduce', type=int, default=5, help='factor by which the magnetogram is reduced')
parser.add_argument('--do_download', type=int, default=0, help='download the files')
args = parser.parse_args()

# get the magnetogram files

big_path, small_path = get_magnetogram_file(cr=args.cr, datdir=args.datdir, force_download=args.do_download, reduce=args.reduce)

write_params_file(args.datdir, args.cr, small_path, args.reduce)

