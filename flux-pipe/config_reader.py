"""
configReader.py - A module for loading configurations from a .ini file and setting environment variables accordingly.
=======================================================================================================================

The main function is `get_configs`, which reads a configuration file,
sets environment variables according to the configuration, and optionally prints the configuration.

The `delete_lines` function is used to delete lines from a file that contain any of a list of strings.

The `set_env` function checks and sets an environment variable.

The `print_config` and `print_one_config` functions are helper functions to print configurations.

Example
-------
To use this module, you can import it and call the `get_configs` function:

    import config_reader
    configs = config_reader.get_configs("DEFAULT")

This will load the configuration named "DEFAULT" from the default .ini file,
set the environment variables, and return the configuration as a dictionary.
"""

import configparser
import os, sys
import ast
import json

global init
init = False

# Initialize default directories and filenames
default_dir = os.path.join(os.environ.get('FL_PREFIX'), "flux-pipe/config")
vars_filename = os.path.join(default_dir, "variables.json")
config_filename = os.path.join(default_dir, "config.ini")
selected_section = "DEFAULT2"



#################################################################
# Config Functions

def print_config(config, config_name=None):
    """
    Print a configuration.

    Parameters
    ----------
    config : ConfigParser
        The loaded configurations
    config_name : str, optional
        The name of the configuration to print, by default None
    """
    # If no specific configuration is specified, print all configurations
    if config_name is None:
        print("\nAll configurations:\n~~~~~~~~~~~~~~~~~~~~~~~")
        for section in config.items():
            sect = section[1]
            print(section[0])
            for key in sect:
                print(f"\t{key}: {sect[key]}")
            print("")
    else:
        print_one_config(config)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

def print_one_config(config):
    """
    Print a specific configuration.

    Parameters
    ----------
    config : dict
        The loaded configuration as a dictionary
    """
    print('\n(py) Configuration loaded!\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

    name = config['config_name']
    print(f'\n[{name}]')
    for key in sorted(config):
        print(f"\t{key}\t  =  {config[key]}")
    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

def delete_lines(original_file, strings_to_delete, temp_file='temp.txt'):
    """
    Delete lines from a file that contain any of a list of strings.

    Parameters
    ----------
    original_file : str
        The path to the file
    strings_to_delete : list
        The strings to delete lines for
    temp_file : str, optional
        The path to a temporary file, by default 'temp.txt'
    """
    if not isinstance(strings_to_delete, list):
        strings_to_delete = [strings_to_delete]

    with open(original_file, 'r', encoding="utf-8") as file, open(temp_file, 'w', encoding="utf-8") as temp:
        for line in file:
            if not any(string in line for string in strings_to_delete):
                temp.write(line)

    os.remove(original_file)
    os.rename(temp_file, original_file)


def get_configs(config_name=None, config_file='flux-pipe/config/config.ini', print_one=False, verbose=False, silent=False):
    """
    Load configurations from a .ini file and set environment variables accordingly.

    Parameters
    ----------
    config_name : str, optional
        The name of the configuration to load, by default "DEFAULT"
    config_file : str, optional
        The path to the .ini file, by default 'flux-pipe/config/config.ini'
    print_one : bool, optional
        Whether to print the loaded configuration, by default True
    verbose : bool, optional
        Whether to print verbose output, by default False

    Returns
    -------
    dict
        The loaded configuration as a dictionary
    """
    if silent:
        verbose = False

    initialize_files()

    #################################################################
    # Begin Loading Configs
    # print("verbose = ", verbose)

    if verbose: print("\n\nReading Configuration File...", end="")
    configs = configparser.ConfigParser()

    fl_prefix = os.environ["FL_PREFIX"]
    abs_config_file = os.path.join(fl_prefix, config_file)

    file_exists = os.path.exists(abs_config_file)
    if not file_exists:
        print(f"\n\tNo config file found at {abs_config_file}\n")
        sys.exit(1)

    # print("verbose =  ", verbose)

    # if verbose:
    #     # print("\n")
    #     print("\n\tFile Path: ", abs_config_file)
    #     print("\n\tFile exists: ", file_exists)

    configs.read(abs_config_file)

    if len(configs.sections()) < 1:
        print(f"\n\tNo sections found in {abs_config_file}\n")
        sys.exit(1)

    if verbose:
        print("\n\tAvailable sections:", configs.sections())

    try:
        if config_name is None:
            config_name = dict(configs["DEFAULT"])['config_name']
        if not silent:
            if not (verbose or print_one):
                print(f"\n\nLoading Configuration [{config_name}]...", end="")
        the_config = dict(configs[config_name])
        # the_config['config_name'] = config_name
    except KeyError:
        print(f"Section '{config_name}' not found in {abs_config_file}. Available sections are {configs.sections()}.")
        sys.exit(1)

    # if verbose:
    #     print("\n\tDone!")

    # Do some string formatting and calculations
    the_config['abs_rc_path'] = RC_FILE_PATH = os.path.expanduser(the_config['rc_path'])
    the_config["run_script"]    = os.path.join(the_config['fl_prefix'], the_config["run_script"])
    the_config["rotations"]     = rots = ast.literal_eval(the_config["rotations"])
    the_config["fluxon_count"]  = nflux= ast.literal_eval(the_config["fluxon_count"])
    the_config["n_jobs"]        = str(len(rots) * len(nflux))

    # Set the perl and flux paths
    set_env('RC_PATH',    the_config['rc_path'],      RC_FILE_PATH)
    set_env('PL_PREFIX',  the_config['pl_prefix'],    RC_FILE_PATH)
    set_env('FL_PREFIX',  the_config['fl_prefix'],    RC_FILE_PATH)
    set_env('WHICH_CONFIG',  config_name,    RC_FILE_PATH)


    if not silent:
        if (verbose or print_one):
            print_one_config(the_config)
        else:
            print(" Done!\n\n")

    return dict(the_config)

# Function to select config file
def select_config_file(filename, section):
    global config
    config = configparser.ConfigParser()
    config.read(filename)
    update_section(section)

# Function to update the current section
def update_section(new_section):
    global selected_section
    if new_section in config.sections():
        selected_section = new_section
    else:
        print(f"Section not found: {new_section}")




#################################################################
# Environment Functions

def set_env(key, value, RC_FILE_PATH, verbose=False):
    """
    Check and set an environment variable.

    Parameters
    ----------
    key : str
        The name of the environment variable
    value : str
        The value of the environment variable
    RC_FILE_PATH : str
        The path to the RC file
    verbose : bool, optional
        Whether to print verbose output, by default False
    """
    if key not in os.environ or os.environ[key] != value:
        delete_lines(RC_FILE_PATH, key)
        os.environ[key] = value
        with open(RC_FILE_PATH, 'a+', encoding="utf-8") as file:
            file.write(f"export {key}={value}\n")
            if verbose:
                print(f"\tSet {key} to {os.environ[key]}")


def get_env(key, default=None):
    """
    Get an environment variable.

    Parameters
    ----------
    key : str
        The name of the environment variable
    default : str, optional
        The default value of the environment variable, by default None

    Returns
    -------
    str
        The value of the environment variable
    """
    return os.environ.get(key, default)

def get_envs(verbose=False, silent=False):
    """
    List all environment variables and their values.
    """

    initialize_files()
    if not silent: print("\nLoading Envs... Done!")
    if verbose and not silent:
        print("\nEnvironment Variables:\n~~~~~~~~~~~~~~~~~~~~~~~\n")
        print_all_envs()

    return os.environ.copy()

def print_all_envs():
    # Find the length of the longest key for alignment
    max_key_length = max(len(key) for key in os.environ.keys())

    for key, value in sorted(os.environ.items()):
        if ":" in value:
            vals = value.split(":")
            print(f"\n{key.rjust(max_key_length)}  =  {vals.pop(0)}")
            for val in vals:
                if val:
                    print(max_key_length*" ", "  >", val)
            print("")
        else:
            print(f"{key.rjust(max_key_length)}  =  {value}")




#################################################################
# Variables Functions


# Function to select variables file
def select_vars_file(filename):
    global vars_filename
    vars_filename = filename
    ensure_json_file_exists()

# Function to ensure JSON file exists
def ensure_json_file_exists():
    if not os.path.exists(vars_filename):
        with open(vars_filename, 'w') as f:
            json.dump({}, f)

# Helper function to read JSON data from a file
def read_json_file(filename):
    try:
        with open(filename, 'r') as f:
            json_text = f.read()
        return json.loads(json_text)
    except FileNotFoundError:
        raise Exception(f"Could not open file '{filename}'")

# Helper function to write JSON data to a file
def write_json_file(data, filename=None):
    if filename is None:
        filename = vars_filename
    try:
        with open(filename, 'w') as f:
            f.write(json.dumps(data))
    except FileNotFoundError:
        raise Exception(f"Could not open file '{filename}'")

# Helper function to ensure JSON file exists
def ensure_json_file_exists():
    if not os.path.exists(vars_filename):
        write_json_file({})

# Function to select variables file
def select_vars_file(filename):
    global vars_filename
    vars_filename = filename
    ensure_json_file_exists()

# Function to delete all variables from JSON file
def delete_all_vars(should_print=False):
    write_json_file({})
    if should_print:
        print("\nAll variables have been deleted from the JSON file.\n")

# Function to set a variable
def set_variable(var_name, new_value):
    ensure_json_file_exists()
    data = read_json_file(vars_filename)
    data[var_name] = new_value
    write_json_file(data)
    return new_value

# Function to get a variable
def get_variable(var_name):
    ensure_json_file_exists()
    data = read_json_file(vars_filename)
    if var_name not in data:
        raise Exception(f"Variable not found: {var_name}")
    return data[var_name]

# Function to print all variables
def print_all_variables(data):
    print("\nAll variables from JSON file:\n")
    if data:
        max_key_length = max(len(key) for key in data.keys())
        for key in sorted(data.keys()):
            print(f"\t{key.ljust(max_key_length)} : {data[key]}")
    else:
        print("\tNo variables found.\n")
    print()

# Function to get variable values from JSON file
def get_vars(verbose=False, silent=False):
    initialize_files()

    if not silent:
        print("\nLoading Vars...", end="")
    ensure_json_file_exists()
    data = read_json_file(vars_filename)
    if verbose and not silent:  # and configs['verbose']:
        print_all_variables(data)
    # if not silent:
    #     print(" Done!\n")
    return data


#################################################################
# Helper Functions

# Function to initialize default files
def initialize_files():
    global init
    if init:
        return
    select_vars_file(vars_filename)

    if os.path.exists(config_filename):
        # First, load the DEFAULT section to get CONFIG_NAME
        select_config_file(config_filename, 'DEFAULT2')
        config_name = config.get('DEFAULT2', 'config_name', fallback='DEFAULT2')
        # Now, load the section named by CONFIG_NAME
        update_section(config_name)
    else:
        print("Please select the config file using 'select_config_file'.")
    init = True

def get_all(verbose=False, silent=False):
    envs    = get_envs(verbose=verbose, silent=silent)
    configs = get_configs(verbose=verbose, silent=silent)
    varbs   = get_vars(verbose=verbose, silent=silent)
    return configs, varbs, envs

def test_reader(verbose = True, silent=False):
    delete_all_vars()
    set_variable("test_var", "test_value")
    configs, varbs, envs = get_all(verbose=verbose, silent=silent)
    return configs, varbs, envs

if __name__ == "__main__":
    # Usage Example
    configs, varbs, envs = get_all(verbose=True, silent=False)
