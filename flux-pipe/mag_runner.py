import subprocess
from tqdm import tqdm
import os

# Carrington Rotations to Treat
# rotations = [2101, 2134, 2135, 2150, 2159, 2192, 2193, 2219]
# rotations = [
#              2101,
#              2110,
#             #  2120,
#              2135,    # very little polar field
#             #   2140,
#             #  2150, 
#              2160,      # mixed polarities at the equator
#             #  2170,
#              2183,    # mixed polarities at the equator
#              2193,    # a good dipole, or polar field
#             #  2200,
#              2210, 
#             #  2215,
#             #  2219, 

#              2231,    # literally only 8 open field lines
#             #  2240,
#              2250, 
#              ]
# rotations = [2130, 2160, 2193, 2219, 2231]
rotations = [2193] #, 2160, 2193, 2219, 2231]
recompute = 0
nflux = 1000
reduction = 2
batch_name = "paperfigs_6"

do_survey = True

# Path to the PDL script
flux_pipe_path = "/Users/cgilbert/vscode/fluxon-mhd/pdl/flux-pipe/"
pdl_script_path = flux_pipe_path + "magnetogram2wind.pdl"
os.chdir(flux_pipe_path)

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
        print("\n    >>verbose = True. All stdout from the processes will be printed following each iteration.<<\n")
else:
    print("\n\nProcessing the following CR: ", rotations, "\n\n")
    print(f"\n\n Batch Name fluxon_{batch_name} \n\n")
    

to_break = 0
# Run the PDL script using subprocess
for rot in tqdm(rotations, desc="    Processing Rotations. Overall progress", unit="rotation"):
    try:
        if do_survey:
            recompute = 0
            
            for reduction in [2]:
                to_break = 0
                if reduction == 2:
                    do_flux = [
                    # 250, 250, 500, 1000, 1500, 2000, 2500, 
                    # 3000, 4000, 5000, 6000, 8000, 10000,
                    # 500, 1000, 2000,

                    1500, 2000, 2500,
                    # 6000, 8000, 10000
                    ]
                else:
                    do_flux = [
                    # 250, 500, 1000, 1500, 2000, 2500, 
                    # 500, 1000, 2000,
                    # 3000, 
                    5000, 
                    6000, 8000, 10000,
                    12000, 14000, 16000
                    ]

                for nflux in do_flux:
                    result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download), str(recompute), str(nflux), str(batch_name)], capture_output=capture)
                    if capture and verbose:
                        print(result.stdout.decode())
                    if result.returncode != 0:
                        to_break += 1
                        if to_break > 2:
                            break
        else:
            result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download), str(recompute), str(nflux), str(batch_name)], capture_output=capture)
            if capture and verbose:
                print(result.stdout.decode())
            if result.returncode != 0:
                to_break += 1
                if to_break > 2:
                    break
    except Exception as e:
        print(e)
    



