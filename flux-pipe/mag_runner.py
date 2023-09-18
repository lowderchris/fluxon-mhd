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
    do_flux:        List of fluxon count values to analyze.
    do_survey:      Flag to run the fluxon analysis on a set of fluxon numbers and/or rotations.
    plot_only:      Flag to skip everything except the wind plotting at the end.
    recompute:      Flag to reperform the fluxon analysis from scratch.
    reduction:      The factor by which the magnetogram is reduced.
    capture:        Flag to capture the stdout of the subprocess.
    verbose:        Flag to print detailed output.
    do_download:    Flag to download the files.

Functions:
    add_top_level_dirs_to_path:  Utility function imported from 'pipe_helper'.

Dependencies:
    subprocess, os, tqdm, pipe_helper

Example:
    Modify the configuration variables as needed and run `python mag_runner.py`.

Author:
    Gilly <gilly@swri.org> (and others!)

Note:
    Each iteration takes around a minute. Please be patient if verbose is set to True.
"""

import subprocess
from os import chdir
import os
from tqdm import tqdm

flux_dir = os.environ.get("FL_PREFIX")
chdir(flux_dir)

pdl_run_script_path = f"{flux_dir}/flux-pipe/magnetogram2wind.pdl"
print(pdl_run_script_path, "\n\n\n\n")


batch_name = "new_test"
rotations = [2108]  #  [2160, 2193, 2219, 2231]
do_flux = [1000]  # , 2000, 3000, 4000, 5000, 6000, 8000, 10000]

# Options
plot_only = 0  # skip everything except the wind plotting at the end
recompute = 0  # reperform the fluxon analysis from scratch
reduction = 2
capture = False
verbose = True
do_download = 0

to_break = 0
print("")
n_jobs = len(rotations) * len(do_flux)

print_run_command = False

with tqdm(total=n_jobs, unit="runs") as pbar:
    # Iterate over the Carrington Rotations
    for rot in rotations:
        # Update the description with the current iteration
        pbar.set_description(f"Rotation {rot}, n_Fluxon {do_flux}")
        print("\n\n\n")

        for nflux in do_flux:
            args = [str(rot), str(reduction), str(do_download), str(recompute), str(nflux), batch_name, str(plot_only)]

            if print_run_command:
                print("Run Command: \n   ", "perl", pdl_run_script_path, " ".join(args), "\n")

            result = subprocess.run(["perl", pdl_run_script_path] + args, check=False)

        # Increment the progress bar
        pbar.update(1)
