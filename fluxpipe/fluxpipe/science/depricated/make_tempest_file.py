import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from fluxpipe.helpers.pipe_helper import configurations
import pandas as pd
# configs = configurations(debug=False)


def parse_args():
    # Create the argument parser
    print("\n\tMaking Tempest File in Python.")
    configs = configurations()
    import argparse
    parser = argparse.ArgumentParser(description=
                            'This script plots the expansion factor of the given radial_bmag_all.dat')
    parser.add_argument('--cr',     type=int, default=configs['rotations'][0],    help='Carrington Rotation')
    parser.add_argument('--dat_dir',type=str, default=configs["data_dir"],        help='data directory')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"],      help='select the batch name')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='magnetogram file')
    parser.add_argument('--show',   type=int, default=0)
    parser.add_argument('--file',   type=str, default=None)
    parser.add_argument('--adapt',  type=int, default=configs["adapt"],           help='Use ADAPT magnetograms')
    args = parser.parse_args()
    filename = args.file or f'{args.dat_dir}/batches/{args.batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    configs = configurations(args=args)
    return filename, configs



# Main Code
if __name__ == "__main__":
    filename, configs = parse_args()

    print(filename)

    #load the table from disk in pandas
    df = pd.read_csv(filename, delim_whitespace=True)

    # print(df["radius", "b_mag"])
    print(df.columns)
    df.groupby("fnum").plot()