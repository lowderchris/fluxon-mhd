"""_summary_

Returns:
    _type_: _description_
"""

import os
import sys
import matplotlib.pyplot as plt

def get_ax(ax=None):
    """ Get the fig and ax, either a new one or the one passed in.

    Args:
        ax (_type_, optional): _description_. Defaults to None.

    Returns:
        fig: pyplot Figure
        ax0: pyplot Axis
    """
    if ax is not None:
        fig, ax0 = ax.get_figure(), ax
    else:
        fig, ax0 = plt.subplots(1)
    return fig, ax0

def add_parent_dir(): #depricated and wrong
    """_summary_
    """
    # Get the parent directory path
    this_cwd = os.getcwd()
    if "plotting" in this_cwd:
        pipdir = os.path.abspath(os.path.join(this_cwd, ".."))
    elif "pipe" in this_cwd:
        pipdir = this_cwd
    else:
        pipdir = os.path.abspath(os.path.join(this_cwd, "flux-pipe"))

    # Add the parent directory to the module search path
    sys.path.append(pipdir)

def add_pipedir(do_print=False):
    """_summary_

    Args:
        do_print (bool, optional): _description_. Defaults to False.
    """
    # get the pipedir environment variable
    pipedir = os.environ.get("PIPEPATH")
    #add the directory to the search path
    sys.path.append(pipedir)
    if do_print:
        print(f"\nAdded {pipedir} to the module search path!\n")


add_pipedir()
