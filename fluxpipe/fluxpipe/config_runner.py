"""
Config Runner: Automate the FLUXpipe Workflow for Multiple Carrington Rotations
===============================================================================

This module serves as a workflow automation tool for the FLUXpipe project. It runs
the `magnetogram2wind.pdl` script for various Carrington Rotations (CRs) specified
in an external 'config.ini' file. The script leverages Python's subprocess module to
invoke the PDL script and uses the tqdm library for tracking progress.

Attributes
----------
configs : dict
    The configuration parameters loaded from 'config.ini'.

Functions
---------
None

Examples
--------
1. Edit the 'config.ini' file to configure your desired settings.
2. Run the script using `python config_runner.py`.

Authors
-------
Gilly <gilly@swri.org> (and others!)

"""

import subprocess
from tqdm import tqdm
from fluxpipe.fluxpipe.pipe_helper import configurations

configs = configurations(debug=False)
def run():
    # Initialize a progress bar with the total number of jobs to run
    with tqdm(total=int(configs["n_jobs"]), unit="runs") as pbar:

        # Loop through specified adapt values
        for adapt in configs["adapts"]:
            configs = configurations(adapt=adapt)

            # Loop through specified Carrington Rotations
            for rot in configs["rotations"]:

                # Loop through specified fluxon counts
                for nflux in configs["fluxon_count"]:

                    # Update the progress bar description
                    pbar.set_description(f"Rotation {rot}, n_fluxon {nflux}")

                    # Execute the PDL script with the current parameters
                    result = subprocess.run(["perl", configs["run_script"],
                                            str(rot), str(nflux), str(adapt)], check=False)

                    # Update the progress bar
                    pbar.update(1)
