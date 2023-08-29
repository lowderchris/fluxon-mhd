"""
Plotting Scalability of Fluxon Run Time
=======================================

This script generates a plot that shows how the run time scales with the number of fluxons.
It provides insights into the performance and scalability of the fluxon simulation for different
configurations and resolutions.

Usage:
    python plot_scalability_single.py

Parameters:
    None (script-based, parameters are hardcoded)

Functions:
    None (script-based)

Example:
    python plot_scalability_single.py

Dependencies:
    os.path, numpy, matplotlib.pyplot, pandas, py_plot_helper

Output:
    A plot saved as both PNG and PDF formats in the specified directory.

Author:
    Gilly <gilly@swri.org> (and others!)

"""


import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from py_plot_helper import load_data

run_name = "fluxon r3"
datdir = "/Users/cgilbert/vscode/fluxons/fluxon-data/"
file_path = os.path.join(datdir, f"{run_name}/pipe_log.txt")


## Real Code

# Load the Array
df = load_data(file_path)
df = df.drop(df.index[0])
print(df)

cr = df.loc['cr'].unique()[0]

n_reductions = len(df.loc["r"].unique())
ny_plots = min(2, n_reductions)
nx_plots = min(int(np.ceil(n_reductions/ny_plots)), n_reductions)
fig, axarray = plt.subplots(nx_plots, ny_plots, figsize=(12, 8), sharex="all", sharey="all")
if not isinstance(axarray, np.ndarray):
    axarray = np.array(axarray)
flaxarray = axarray.flatten()


if len(df.loc["r"].unique()) == 1:
    reduce = int(df.loc["r"].unique()[0])
else:
    reduce = 0

to_plot = ["n_want", "TrOpen", "n_out", 'ttime', "Success"]
independent = "n_actual"
ax=None
for (ax, reduction, rx, ry) in zip(flaxarray, df.loc["r"].unique(),
                                df.loc['rx'].unique(), df.loc['ry'].unique()):
    ax.set_title(f"CR{int(cr)}, Resolution: ({int(rx)}, {int(ry)})")
    ax.set_yscale('log')
    ax.set_xscale('log')


    ax.set_xlabel("Actual Number of Fluxons")
    ax.set_ylabel("Values", labelpad=-5)
    for ii, line_name in enumerate(to_plot):
        this_line = df.loc[df.loc["r"] == reduction][line_name]
        abscissa  = df.loc[df.loc["r"] == reduction][independent]

        # Filter NaN values
        valid_indices = np.logical_and(~np.isnan(abscissa), ~np.isnan(this_line))
        abscissa_valid = abscissa[valid_indices]
        this_line_valid = this_line[valid_indices]
        if abscissa_valid.any():
            coefficients = np.polyfit(abscissa_valid, this_line_valid, 1)
            fit_line = np.poly1d(coefficients)(abscissa)
            round_coeff = [round(ccc, 2) for ccc in coefficients]
            this_label = f"{line_name}: {round_coeff[0]}"
            color = f"C{ii}"

            if line_name in ["Success"]:
                this_line_valid[this_line_valid <= 0.1] = 0.1
                # the same plot but stair-stepped
                ax.step(abscissa_valid, this_line_valid, color=color,
                        where='mid', marker='o', label="Success")
                # find the index of the first non-"1" value in this_line_valid
                non_one_index = next((i for i, val in enumerate(this_line_valid) if val != 1), None)
                non_one_loc = abscissa_valid[non_one_index]
                ax.axvline(non_one_loc, color="k", ls=":")


            elif line_name in ["ttime"]:
                ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label, marker='o')
            else:
                ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label, marker='o')
                this_line_valid[this_line_valid <= 0.1] = 0.1
                fit_line[fit_line <= 0.1] = 0.1
                # fit_line[fit_line == 0] = 0.1
                ax.plot(abscissa, fit_line,  color=color, ls=":")
    ax.axhline(100, color="k", ls="--", label="100 Fluxons")
    ax.legend(loc=(0.6, 0.25 ))


ax.set_ylim((5*10**(-2), 2.1*(10**4)))
ax.set_xlim((80, 15000))

# line that changes the ylim to be from 1 to the current ymax
fig.set_size_inches(5, 5.5)
# Fit line
plt.tight_layout()

img_path = os.path.join(datdir, f"{run_name}/imgs/scalability_{reduce}.png")
img_path_pdf = os.path.join(datdir, f"{run_name}/imgs/scalability_{reduce}.pdf")
plt.tight_layout()

plt.savefig(img_path)
plt.savefig(img_path_pdf)
# plt.show()
plt.close(fig)
