"""
Manually Run FLUXpipe as a script in a loop
=============================================

This module is used to automate the execution of the `magnetogram2wind.pdl` script.
It iterates through a list of Carrington Rotations and fluxon counts, running the script
for each combination. This is manual in the sense that you must specify the parameters directly,
rather than putting them in a configuration file.

Attributes
----------
flux_dir : str
    The directory where the FLUXpipe project is located.
pdl_run_script_path : str
    The path to the `magnetogram2wind.pdl` script.
batch_name : str
    The name of the batch for the current run.
rotations : list of int
    List of Carrington Rotations to process.
do_flux : list of int
    List of fluxon counts to process.

Options
-------
plot_only : int
    Skip everything except the wind plotting at the end.
recompute : int
    Perform the fluxon analysis from scratch.
reduction : int
    The reduction factor for the magnetogram.
capture : bool
    Whether to capture the output of the subprocess.
verbose : bool
    Whether to print verbose output.
do_download : int
    Whether to download new data.

Example
--------
python manual_runner.py

Authors
-------
Gilly <gilly@swri.org> (and others!)

"""

import subprocess
from os import chdir
import os
from tqdm import tqdm

# Initialize variables and options
flux_dir = os.environ.get("FL_PREFIX")
chdir(flux_dir)
pdl_run_script_path = f"{flux_dir}/fluxpipe/magnetogram2wind.pdl"


batch_name = "new_test"
rotations = [2108]
do_flux = [1000]
plot_only = 0
recompute = 0
reduction = 2
capture = False
verbose = True
do_download = 0
print_run_command = False

n_jobs = len(rotations) * len(do_flux)


###############################################################################

# Initialize a progress bar with the total number of jobs to run
with tqdm(total=n_jobs, unit="runs") as pbar:
    # Loop through specified Carrington Rotations
    for rot in rotations:
        # Update the progress bar description
        pbar.set_description(f"Rotation {rot}, n_Fluxon {do_flux}")

        # Loop through specified fluxon counts
        for nflux in do_flux:
            args = [str(rot), str(reduction), str(do_download), str(recompute), str(nflux), batch_name, str(plot_only)]

            # Print the run command if verbose mode is enabled
            if print_run_command:
                print("Run Command: \n   ", "perl", pdl_run_script_path, " ".join(args), "\n")

            # Execute the PDL script with the current parameters
            result = subprocess.run(["perl", pdl_run_script_path] + args, check=False)

        # Update the progress bar
        pbar.update(1)
