import subprocess
from tqdm import tqdm
from fluxpipe.helpers.pipe_helper import configurations
import timeout_decorator
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

configs = configurations(debug=False)

def ask_for_cores():
    max_cores = multiprocessing.cpu_count()
    print(f"Your machine has {max_cores} cores available.")
    try:
        n_cores = int(input("How many cores would you like to use? "))
        if 1 <= n_cores <= max_cores:
            return n_cores
        else:
            print(f"Please enter a number between 1 and {max_cores}.")
            return ask_for_cores()
    except ValueError:
        print("Please enter a valid number.")
        return ask_for_cores()

@timeout_decorator.timeout(600)  # Set a 600-second timeout for each subprocess call
def run_pdl_script(params):
    """
    Executes the PDL script with a timeout.

    Parameters:
    - params: A tuple containing (rot, nflux, adapt)
    """
    rot, nflux, adapt = params
    try:
        subprocess.run(["perl", configs["run_script"], str(rot), str(nflux), str(adapt)], check=False)
        return ("Success", rot, nflux, f"[{rot}, {nflux}, {adapt}]")
    except subprocess.TimeoutExpired:
        return ("Timeout", rot, nflux, f"\033[93m[Timeout {rot}, {nflux}, {adapt}]\033[0m")
    except Exception as e:
        return ("Error", rot, nflux, f"\033[91m[Error {rot}, {nflux}, {adapt}: {e}]\033[0m")

def run_parallel(n_cores):
    tasks = []

    # Prepare tasks
    for adapt in configs["adapts"]:
        for rot in configs["rotations"]:
            for nflux in configs["fluxon_count"]:
                tasks.append((rot, nflux, adapt))

    # Execute tasks in parallel
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = {executor.submit(run_pdl_script, task): task for task in tasks}
        for future in tqdm(as_completed(futures), total=len(tasks), unit="run"):
            result = future.result()
            print(result[3])  # Print the colored and identified result

if __name__ == "__main__":
    n_cores = ask_for_cores()
    run_parallel(n_cores)
