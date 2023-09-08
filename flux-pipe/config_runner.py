"""
Mag Runner: Automate the FLUXpipe Workflow for Multiple Carrington Rotations
============================================================================

This script automates the entire FLUXpipe workflow by executing the
`magnetogram2wind.pdl` script for a set of specified Carrington Rotations (CRs).
It provides options to control various parameters such as the number of fluxons,
the reduction factor, and more. The script utilizes the subprocess module to run
the PDL script and tqdm for progress tracking.

Usage:
    python mag_runner.py [config_name]

Configuration Variables:
    which_config:   The name of the configuration to load from 'config_reader'.
    configs:        The loaded configuration as a dictionary.

Parameters Controlled by Configuration:
    fl_prefix:      The base directory where the FLUXpipe workflow is located.
    rotations:      List of Carrington Rotations to process.
    fluxon_count:   List of fluxon count values to analyze.
    magnetogram_reduction: The factor by which the magnetogram is reduced.
    do_download:    Flag to download the files.
    recompute:      Flag to reperform the fluxon analysis from scratch.
    batch_name:     The name of the batch for the operation.
    plot_only:      Flag to skip everything except the wind plotting at the end.
    n_jobs:         The total number of jobs to be processed.

Dependencies:
    subprocess, os, tqdm, config_reader, sys

Example:
    Run `python mag_runner.py` with an optional configuration name as an argument.

Author:
    Gilly <gilly@swri.org> (and others!)

Note:
    Each iteration takes around a minute. Please be patient.
"""

import subprocess
from os import chdir
from tqdm import tqdm
from config_reader import get_all

# Print a separator for better readability
print("\n|\n|\n|\n|\n|\n|----------------------------------------------------|\n|\n|\n|\n|")

# Load the configuration
configs, varbs, envs = get_all(verbose=True, silent=True)


# Change the current directory to the FLUXpipe directory
chdir(configs["fl_prefix"])

# Initialize the progress bar
with tqdm(total=int(configs["n_jobs"]), unit="runs") as pbar:
    for rot in configs["rotations"]:
        # Iterate over the Carrington Rotations

        for nflux in configs["fluxon_count"]:
            # Iterate over the fluxon counts

            # Update the description with the current iteration
            pbar.set_description(f"Rotation {rot}, n_fluxon {nflux}")

            # Set the run arguments
            args = [str(rot), str(configs["mag_reduct"]), str(configs["do_download"]),
                    str(configs["recompute"]), str(nflux), configs["batch_name"], str(configs["plot_only"])]

            # Execute the perl script
            result = subprocess.run(["perl", configs["run_script"]] + args, check=False)

            # Update the progress bar
            pbar.update(1)
