"""
configReader.py - A module for loading configurations from a .ini file and setting environment variables accordingly.
=======================================================================================================================

The main function is `load_configs`, which reads a configuration file,
sets environment variables according to the configuration, and optionally prints the configuration.

The `delete_lines` function is used to delete lines from a file that contain any of a list of strings.

The `check_and_set` function checks and sets an environment variable.

The `get_rc_path` function gets the path to the RC file from a configuration.

Example
-------
To use this module, you can import it and call the `load_configs` function:

    import config_loader
    configs = config_loader.load_configs("Section1")

This will load the configuration named "Section1" from the default .ini file,
set the environment variables, and return the configuration as a dictionary.
"""

import configparser
import os
import ast


def load_configs(config_name="DEFAULT", config_file='flux-pipe/config/config.ini', \
                 print_one=True, verbose=False):
    """
    Load configurations from a .ini file and set environment variables accordingly.

    Parameters
    ----------
    config_name : str, optional
        The name of the configuration to load, by default "DEFAULT"
    config_file : str, optional
        The path to the .ini file, by default 'fluxon-mhd/flux-pipe/config/config.ini'
    print_one : bool, optional
        Whether to print the loaded configuration, by default True
    print_all : bool, optional
        Whether to print all configurations, by default False
    verbose : bool, optional
        Whether to print verbose output, by default False

    Returns
    -------
    dict
        The loaded configuration as a dictionary
    """

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
            print_one_config(config, config_name)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    def print_one_config(config):
        """
        Print a configuration.

        Parameters
        ----------
        config : ConfigParser
            The loaded configurations
        config_name : str, optional
            The name of the configuration to print, by default None
        """
        name = config['config_name']
        print(f'\n>>{name}<< Configuration loaded!\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        for key in config:
            print(f"\t{key} :   {config[key]}")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


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

    def check_and_set(key, value, RC_FILE_PATH, verbose=False):
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
        """
        if key not in os.environ or os.environ[key] != value:
            delete_lines(RC_FILE_PATH, key)
            os.environ[key] = value
            with open(RC_FILE_PATH, 'a+', encoding="utf-8") as file:
                file.write(f"export {key}={value}\n")
                if verbose:
                    print(f"\tSet {key} to {os.environ[key]}")

    #################################################################
    # End Helper Functions

    #################################################################
    # Begin Loading Configs

    if verbose: print("\n\nReading Configuration File...", end="")
    configs = configparser.ConfigParser()
    # /Users/cgilbert/vscode/fluxons/fluxon-mhd
    print()
    fl_prefix = os.environ["FL_PREFIX"]
    abs_config_file = os.path.join(fl_prefix, config_file)

# /Users/cgilbert/vscode/fluxons/fluxon-mhd/flux-pipe/config/config.ini
# /Users/cgilbert/vscode/fluxons/fluxon-mhd/flux-pipe/fluxon-mhd/flux-pipe/config/config.ini


    file_exists = os.path.exists(abs_config_file)
    if not file_exists:
        print(f"\n\tNo config file found at {abs_config_file}\n")
        exit(1)

    print("\n\n")
    print("Reading config file:", abs_config_file)
    print("")
    print("Config file exists:", file_exists)


    configs.read(abs_config_file)

    if len(configs.sections()) < 1:
        print(f"\n\tNo sections found in {abs_config_file}\n")
        exit(1)

    print("Available sections:", configs.sections(), "\n")

    try:
        the_config = dict(configs[config_name])
    except KeyError:
        print(f"Section '{config_name}' not found in {abs_config_file}. Available sections are {configs.sections()}.")
        exit(1)

    # the_config = dict(configs[config_name])
    the_config['config_name'] = config_name

    if verbose:
        print("Done!")
        print_config(the_config)
        print('\n\nSetting perl and flux paths...')


    the_config['abs_rc_path'] = RC_FILE_PATH = os.path.expanduser(the_config['rc_path'])
    check_and_set('RC_PATH', the_config['rc_path'], RC_FILE_PATH)
    check_and_set('PL_PREFIX', the_config['pl_prefix'], RC_FILE_PATH)
    check_and_set('FL_PREFIX', the_config['fl_prefix'], RC_FILE_PATH)

    # Do some string formatting and calculations
    the_config["run_script"] = os.path.join(the_config['fl_prefix'], the_config["run_script"])

    def parse_string(s):
        return ast.literal_eval(s)

    rots = the_config["rotations"] = parse_string(the_config["rotations"])
    nflux = the_config["fluxon_count"] = parse_string(the_config["fluxon_count"])
    the_config["n_jobs"] = str(len(rots) * len(nflux))


    #.format()

    if print_one: print_one_config(the_config)

    return dict(the_config)

if __name__ == "__main__":
    # Usage Example
    configs = load_configs("Child1")
