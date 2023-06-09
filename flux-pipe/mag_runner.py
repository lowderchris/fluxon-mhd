import subprocess
from tqdm import tqdm
import os
import sys



rotations = [2193] #, 2160, 2193, 2219, 2231]
recompute = 0
nflux = 1000
reduction = 2
batch_name = "default_batch"

do_survey = False

# Path to the PDL script
flux_pipe_dir = "/Users/cgilbert/vscode/fluxon-mhd/flux-pipe/"
pdl_script_path = flux_pipe_dir + "magnetogram2wind.pdl"
os.chdir(flux_pipe_dir)

# Get the parent directory path
# parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
plot_dir = os.path.abspath(os.path.join(flux_pipe_dir, "plotting"))
# Add the parent directory to the module search path
sys.path.append(plot_dir)


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

            # Increment the progress bar
            pbar.update(1)

        except Exception as e:
            print(e)
        



