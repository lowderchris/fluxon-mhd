import subprocess
from tqdm import tqdm

# Carrington Rotations to Treat
# rotations = [2101, 2134, 2135, 2150, 2159, 2192, 2193, 2219]
rotations = [2183, 2160, 2135, 2130, 2193, 2231] #, 2135, 2160]

# 2183, 2160 have mixed polarities at the equator
# 2135, 2130 have very little polar field
# 2193,      have a good dipole, or polar field
# 2231,      has literally 8 open field lines

# Path to the PDL script
pdl_script_path = "magnetogram2wind.pdl"

capture = False
verbose = True
do_download = 0 # Depricated
reduction = 5
recompute = 1
nflux = 800

if capture:
    print("\n\n\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Magnetogram 2 Wind: Run the entire fluxon pipeline on a set of Carrington rotations.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"    Target Rotations: {rotations}")
    print("    Each iteration takes around a minute. Please be patient.")
    if verbose:
        print("\n    >>verbose = True. All stdout from the processes will be printed following each iteration.<<\n")

# Run the PDL script using subprocess
for rot in tqdm(rotations, desc="    Processing Rotations. Overall progress", unit="rotation"):
    try:
        result = subprocess.run(["perl", pdl_script_path, str(rot), str(reduction), str(do_download), str(recompute), str(nflux)], capture_output=capture)
    except Exception as e:
        print(e)
    if capture and verbose:
        print(result.stdout.decode())



