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
    mhd_dir:        The directory containing the MHD (Magnetohydrodynamics) code.
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
import py_pipe_helper as ph

base_dir = "/Users/cgilbert/vscode/fluxons"
mhd_dir = f"{base_dir}/fluxon-mhd"
flux_pipe_dir = f"{mhd_dir}/flux-pipe"
chdir(mhd_dir)
new_path = ph.add_top_level_dirs_to_path(mhd_dir)


pdl_script_path = f"{flux_pipe_dir}/magnetogram2wind.pdl"

batch_name = "new_test"
# rotations = [2160, 2193, 2219, 2231]
rotations = [2100]
do_flux = [1000]  # , 2000, 3000, 4000, 5000, 6000, 8000, 10000]
do_survey = True  # run the fluxon analysis on a set of fluxon numbers and/or rotations

plot_only = 0  # skip everything except the wind plotting at the end
recompute = 0  # reperform the fluxon analysis from scratch
# nflux = 500
reduction = 2

# Options
capture = False
verbose = True
do_download = 0


if capture:
    print("\n\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Magnetogram 2 Wind: Run the entire fluxon pipeline on a set of Carrington rotations.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"    Target Rotations: {rotations}")
    print("    Each iteration takes around a minute. Please be patient.")
    if verbose:
        print("\n    >>verbose = True. All stdout from the processes \
              will be printed following each iteration.<<\n")
else:
    pass


to_break = 0
print("")
with tqdm(total=len(rotations), unit="rotation") as pbar:
    for rot in rotations:
        try:
            # Update the description with the current iteration
            pbar.set_description(f"Processing Rotation {rot}")
            print("\n")

            if do_survey:
                recompute = 0
                # do_flux = [8000]

                for nflux in do_flux:
                    result = subprocess.run(["perl", pdl_script_path, str(rot),
                            str(reduction), str(do_download), str(recompute),
                            str(nflux), str(batch_name), str(plot_only)],
                            capture_output=capture, check=False)
                    # exit()
                    if capture and verbose:
                        print(result.stdout.decode())
                    if result.returncode != 0:
                        to_break += 1
                        if to_break > 2:
                            break
            else:
                result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction),
                            str(do_download), str(recompute), str(nflux), str(batch_name)],
                            capture_output=capture, check=False)
                if capture and verbose:

                    print(result.stdout.decode())
                if result.returncode != 0:
                    to_break += 1
                    if to_break > 2:
                        break

            # Increment the progress bar
            pbar.update(1)

        except Exception as e:
            print(e)
            assert False, "Error in the main loop"
