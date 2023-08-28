""" This script plots the results of a scalability test.
"""

import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from py_plot_helper import load_data

datdir = "/Users/cgilbert/vscode/fluxons/fluxon-data/"
batch = "scalability_test"
file_path = os.path.join(datdir, "batches", batch, "pipe_log.txt")


## Real Code
# Load the Array
df_scaling_array = load_data(file_path)
df_scaling_array = df_scaling_array.drop(df_scaling_array.index[0])
print(type(df_scaling_array))

n_reductions = len(df_scaling_array.loc["r"].unique())
ny_plots = min(2, n_reductions)
nx_plots = min(int(np.ceil(n_reductions/ny_plots)), n_reductions)
fig, axarray = plt.subplots(nx_plots, ny_plots, figsize=(12, 8), sharex="all", sharey="all")
if not isinstance(axarray, np.ndarray):
    axarray = np.array(axarray)
flaxarray = axarray.flatten()

to_plot = ["n_want", "TrOpen", "n_out", 'ttime', "Success"]
independent = "n_actual"
for (ax, reduction, rx, ry) in zip(flaxarray,
            df_scaling_array.loc["r"].unique(),
            df_scaling_array.loc['rx'].unique(),
            df_scaling_array.loc['ry'].unique()):
    ax.set_title(f"Magnetic Field Map Resolution: ({rx}, {ry})")
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=False))

    ax.set_xlabel("Actual Number of Fluxons")
    ax.set_ylabel("Values")
    for ii, line_name in enumerate(to_plot):
        this_line = df_scaling_array.loc[df_scaling_array.loc["r"] == reduction][line_name]
        abscissa  = df_scaling_array.loc[df_scaling_array.loc["r"] == reduction][independent]

        # Filter NaN values
        valid_indices = np.logical_and(~np.isnan(abscissa), ~np.isnan(this_line))
        abscissa_valid = abscissa[valid_indices]
        this_line_valid = this_line[valid_indices]
        if abscissa_valid.any():
            coefficients = np.polyfit(abscissa_valid, this_line_valid, 1)
            fit_line = np.poly1d(coefficients)(abscissa_valid)
            round_coeff = [round(ccc, 2) for ccc in coefficients]
            this_label = f"{line_name}: {round_coeff[0]}"
            color = f"C{ii}"

            if line_name in ["Success"]:
                this_line_valid[this_line_valid <= 0.1] = 0.1
                # ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label)
                # the same plot but stair-stepped
                ax.step(abscissa_valid, this_line_valid, color=color, marker="o", where='mid')
            elif line_name in ["ttime"]:
                ax.plot(abscissa_valid, this_line_valid, color=color, marker="o", label=this_label)
            else:
                ax.plot(abscissa_valid, this_line_valid, color=color, marker="o", label=this_label)
                this_line_valid[this_line_valid <= 0.1] = 0.1
                fit_line[fit_line <= 0.1] = 0.1
                # fit_line[fit_line == 0] = 0.1
                ax.plot(abscissa_valid, fit_line,  color=color, ls=":", marker="o")
    ax.legend()

# line that changes the ylim to be from 1 to the current ymax
fig.suptitle("Scalability of Fluxon Algorithm")
# Fit line
plt.tight_layout()
img_path = os.path.join(datdir, f"batches/{batch}/imgs/scalability.png")
plt.savefig(img_path)
# plt.show()
plt.close(fig)
