import argparse
from magnetoget import get_magnetogram_files, reduce_fits_image
magneto_file = None

# create the argument parser
parser = argparse.ArgumentParser(description='This script downloads a magnetogram for a particular Carrington Rotation')
parser.add_argument('--cr', type=int, default=2219, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/Fluxon-Scripts-Gilly', help='data directory')
parser.add_argument('--reduce', type=int, default=10, help='factor by which the magnetogram is reduced')
parser.add_argument('--do_download', type=int, default=1, help='download the files')
args = parser.parse_args()

# get the magnetogram files
(hmi_path, mdi_path) = get_magnetogram_files(cr=args.cr, date=None, data_dir=args.dat_dir, do_download=args.do_download)


# reduce the FITS image
if args.reduce:
    magneto_path = reduce_fits_image(hmi_path, target_resolution=None, reduction_amount=args.reduce)
    magneto_file = magneto_path.replace(args.dat_dir, "")
    if magneto_file[0] == "/":
        magneto_file = magneto_file[1:]

# write the parameter file
if magneto_file is not None:
    params_path = args.dat_dir + "/magnetic_target.params"
    with open(params_path, 'w') as fp:
        fp.write("## CR_int, Filename_str, Adapt_bool, Doplot_bool, reduction ##\n")
        fp.write(str(args.cr)+"\n")
        fp.write(str(magneto_file)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(args.reduce))
