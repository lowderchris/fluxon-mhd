import subprocess
from tqdm import tqdm
import os
import sys
import pipe_helper as ph


flux_pipe_dir = "/Users/cgilbert/vscode/fluxon-mhd/fluxpipe/fluxpipe/"
batch_name = "adapt_test2"
# rotations = [2160, 2193, 2219, 2231]
rotations = [2193]
do_flux = [750] #, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
do_survey = True # run the fluxon analysis on a set of fluxon numbers and/or rotations

ADAPT = 1

plot_only = 0 # skip everything except the wind plotting at the end
recompute = 0 # reperform the fluxon analysis from scratch
# nflux = 500
reduction = 2 # if not ADAPT else 'A'
pdl_script_path = ph.add_paths(flux_pipe_dir)

# Options
capture = False
verbose = True
do_download = 0
plot_only = 0

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
    pass
    # print("\n\nProcessing the following CR: ", rotations, "\n\n")
    # print(f"\n\n Batch Name {batch_name} \n\n")


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
                    result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download),
                                             str(recompute), str(nflux), str(batch_name), str(plot_only), str(ADAPT)], capture_output=capture)
                    # exit()
                    if capture and verbose:
                        print(result.stdout.decode())
                    if result.returncode != 0:
                        to_break += 1
                        if to_break > 2:
                            break
            else:
                result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download),
                                         str(recompute), str(nflux), str(batch_name), str(plot_only), str(ADAPT)], capture_output=capture)
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





