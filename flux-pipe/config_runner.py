"""
Mag Runner: Automate the FLUXpipe Workflow for Multiple Carrington Rotations
============================================================================

This script sets the entire FLUXpipe workflow in motion by running the
`magnetogram2wind.pdl` script for a set of specified Carrington Rotations (CRs).
It provides options to control the number of fluxons, the reduction factor, and
other parameters. The script uses the subprocess module to execute the PDL script
and tqdm for progress tracking.

Usage:
    python mag_runner.py

Configuration:
    base_dir:       The base directory where the FLUXpipe workflow is located.
    flux_dir:        The directory containing the MHD (Magnetohydrodynamics) code.
    flux_pipe_dir:  The directory containing the FLUXpipe code.
    batch_name:     The name of the batch for the operation.
    rotations:      List of Carrington Rotations to process.
    fluxon_count:        List of fluxon count values to analyze.
    do_survey:      Flag to run the fluxon analysis on a set of fluxon numbers and/or rotations.
    plot_only:      Flag to skip everything except the wind plotting at the end.
    recompute:      Flag to reperform the fluxon analysis from scratch.
    reduction:      The factor by which the magnetogram is reduced.
    capture:        Flag to capture the stdout of the subprocess.
    verbose:        Flag to print detailed output.
    do_download:    Flag to download the files.

Functions:
    add_top_level_dirs_to_path:  Utility function imported from 'py_pipe_helper'.

Dependencies:
    subprocess, os, tqdm, py_pipe_helper

Example:
    Modify the configuration variables as needed and run `python mag_runner.py`.

Author:
    Gilly <gilly@swri.org> (and others!)

Note:
    Each iteration takes around a minute. Please be patient if verbose is set to True.
"""

import subprocess
from os import chdir
from tqdm import tqdm
import config_reader



which_config = "Child1"



configs = config_reader.load_configs(which_config)
chdir(configs["fl_prefix"])

with tqdm(total=int(configs["n_jobs"]), unit="runs") as pbar:
    for rot in configs["rotations"]:
        # Iterate over the Carrington Rotations

        for nflux in configs["fluxon_count"]:
            # Iterate over the fluxon counts

            # Update the description with the current iteration
            pbar.set_description(f"Rotation {rot}, n_Fluxon {nflux}")
            print("\n\n\n")

            # Set the run arguments
            args = [str(rot), str(configs["magnetogram_reduction"]), str(configs["do_download"]),
                    str(configs["recompute"]), str(nflux), configs["batch_name"], str(configs["plot_only"])]

            # Run the perl script
            result = subprocess.run(["perl", configs["run_script"]] + args, check=False)

            # Increment the progress bar
            pbar.update(1)
