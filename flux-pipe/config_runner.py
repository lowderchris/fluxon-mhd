"""
Config Runner: Automate the FLUXpipe Workflow for Multiple Carrington Rotations
===============================================================================

This script serves as a workflow automation tool for the FLUXpipe project. It runs
the `magnetogram2wind.pdl` script for various Carrington Rotations (CRs) specified
in an external 'config.ini' file. The script leverages Python's subprocess module to
invoke the PDL script and uses the tqdm library for tracking progress.

Usage:
    Run the script using `python mag_runner.py`.

Configuration:
    Edit the 'config.ini' file to set the required configuration parameters.

Dependencies:
    - subprocess
    - os
    - tqdm
    - configparser

Example:
    1. Edit the 'config.ini' to configure your desired settings.
    2. Run the script using `python mag_runner.py`.

Author:
    Gilly <gilly@swri.org> (and others!)

"""

import subprocess
from tqdm import tqdm
import configparser
import os
import ast

def configurations(config_name=None, config_filename="config.ini", debug=False):
    """
    Reads and sanitizes configuration settings from a specified config file.

    Args:
        config_name (str, optional): The specific configuration section to read.
                                     If not provided, defaults to the DEFAULT section.
        config_filename (str, optional): The filename of the configuration file.
                                         Defaults to 'config.ini'.
        debug (bool, optional): Whether to print debug information. Defaults to False.

    Returns:
        dict: Configuration settings as key-value pairs.
    """
    config_obj = configparser.ConfigParser()
    config_path = f"fluxon-mhd/flux-pipe/config/{config_filename}"

    # Search for the configuration file in the current directory and subdirectories
    if not os.path.exists(config_path):
        found = False
        for root, dirs, files in os.walk(os.getcwd()):
            if config_filename in files:
                config_path = os.path.join(root, config_filename)
                found = True
                break
        if not found:
            raise FileNotFoundError("Configuration file not found.")

    # Clean the file content: remove comments and trailing whitespaces
    with open(config_path, "r") as f:
        lines = f.readlines()
    clean_lines = [line.split("#")[0].rstrip() for line in lines]
    clean_content = "\n".join(clean_lines)

    # Parse the clean configuration string
    config_obj.read_string(clean_content)

    # Fallback to section defined in the DEFAULT section if no specific section is provided
    if config_name is None:
        config_name = config_obj["DEFAULT"]['config_name']

    # Extract and further process configuration settings
    the_config = dict(config_obj[config_name])
    the_config['abs_rc_path'] = os.path.expanduser(the_config['rc_path'])
    the_config["run_script"] = os.path.join(the_config['fl_prefix'], the_config["run_script"])
    the_config["rotations"] = ast.literal_eval(the_config["rotations"])
    the_config["fluxon_count"] = ast.literal_eval(the_config["fluxon_count"])
    the_config["n_jobs"] = str(len(the_config["rotations"]) * len(the_config["fluxon_count"]))

    if debug:
        print("\nConfiguration file values:\n--------------------------------")
        for key, value in the_config.items():
            print(f"{key}: \t{value}")
        print("--------------------------------\n\n")

    return the_config


# Load configurations from disk
configs = configurations(debug=True)

# Initialize a progress bar with the total number of jobs to run
with tqdm(total=int(configs["n_jobs"]), unit="runs") as pbar:
    # Loop through specified Carrington Rotations
    for rot in configs["rotations"]:
        # Loop through specified fluxon counts
        for nflux in configs["fluxon_count"]:
            # Update the progress bar description
            pbar.set_description(f"Rotation {rot}, n_fluxon {nflux}")

            # Execute the PDL script with the current parameters
            result = subprocess.run(["perl", configs["run_script"],
                                     str(rot), str(nflux)], check=False)

            # Update the progress bar
            pbar.update(1)
