
import os, sys

def add_parent_dir(): #depricated and wrong
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
    # get the pipedir environment variable
    pipedir = os.environ.get("PIPEPATH")
    #add the directory to the search path
    sys.path.append(pipedir)
    if do_print:
        print(f"\nAdded {pipedir} to the module search path!\n")


add_pipedir()