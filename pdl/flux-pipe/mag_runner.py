import subprocess
from tqdm import tqdm
import os

# Carrington Rotations to Treat
# rotations = [2101, 2134, 2135, 2150, 2159, 2192, 2193, 2219]
rotations = [
             2101,
             2110,
            #  2120,
             2135,    # very little polar field
            #  2140,
            #  2150, 
             2160,      # mixed polarities at the equator
            #  2170,
             2183,    # mixed polarities at the equator
             2193,    # a good dipole, or polar field
            #  2200,
             2210, 
            #  2215,
            #  2219, 
             2231,    # literally only 8 open field lines
            #  2240,
             2250, 
             ]
recompute = 0


# Path to the PDL script
flux_pipe_path = "/Users/cgilbert/vscode/fluxon-mhd/pdl/flux-pipe/"
pdl_script_path = flux_pipe_path + "magnetogram2wind.pdl"
os.chdir(flux_pipe_path)

# Options
nflux = 1500
capture = False
verbose = True
do_download = 0
reduction = 5



if capture:
    print("\n\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Magnetogram 2 Wind: Run the entire fluxon pipeline on a set of Carrington rotations.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"    Target Rotations: {rotations}")
    print("    Each iteration takes around a minute. Please be patient.")
    if verbose:
        print("\n    >>verbose = True. All stdout from the processes will be printed following each iteration.<<\n")
else:
    print("\n\nProcessing the following CR: ", rotations, "\n\n")
    
# Run the PDL script using subprocess
for rot in tqdm(rotations, desc="    Processing Rotations. Overall progress", unit="rotation"):
    try:
        result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download), str(recompute), str(nflux)], capture_output=capture)
    except Exception as e:
        print(e)
        # result = subprocess.run(["perl", os.path.join("/pdl/flux-pipe",pdl_script_path), str(rot), str(reduction), str(do_download), str(recompute), str(nflux)], capture_output=capture)

    if capture and verbose:
        print(result.stdout.decode())
    # exit()



