import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from fluxpipe.helpers.pipe_helper import configurations
configs = configurations()
mpl.use('qt5Agg')
# Specify the path to your CSV file
csv_file_path = 'fluxon-data/seq_fss_theta_mcg2011.csv'

# Specify your column names
column_names = ['Fid', 'fss', 'theta', 'speed']

# Read the CSV file
df = pd.read_csv(csv_file_path, header=None, names=column_names, sep=", ")


vmin=configs["all_vmin"]
vmax=configs["all_vmax"]
# Set the colormap
cmap = plt.cm.autumn #plt.cm.magma #plt.cm.autumn

# Display the DataFrame
# print(df)
fig, [sc_axis, hist_ax] = plt.subplots(2, figsize=(6,6),
                                       gridspec_kw={'height_ratios': [3, 1]})

# Plotting
df.plot(ax=sc_axis, kind='scatter', x='theta', y='fss', c='speed',
        alpha=0.5, edgecolor="none", colormap=cmap, colorbar=False,
        vmin=vmin, vmax=vmax)
sc_axis.set_title(r'Scatter Plot of $f_{ss}$ vs $\theta_b$, ' + str(len(df)) + " points" )
sc_axis.set_xlabel(r'Edge Angle $\theta_b$ [degrees]')
sc_axis.set_xlim((0, 80))
sc_axis.set_ylim((0, 100))
sc_axis.set_ylabel('Expansion Factor $f_{ss}$')
sc_axis.axvline(6, ls="--")
sc_axis.axvline(10)
# plt.show()

# Histogram
def plot_hist(vel_positive, n_bins, vmin, vmax, cmap, hist_ax):
    import numpy as np

    # Calculate the histogram
    hist, bin_edges = np.histogram(vel_positive, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot the histogram
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    for i in range(len(bin_edges) - 1):
        color = cmap(norm(bin_centers[i]))
        hist_ax.bar(
            bin_edges[i], hist[i], width=bin_edges[i + 1] - bin_edges[i], color=color
        )



    # plt.figure(figsize=(10, 6))  # Create a new figure for the histogram
    # hist_ax.hist(df['speed'], bins=20, alpha=0.7, edgecolor='black')
    hist_ax.set_title('Histogram of Speeds')
    hist_ax.set_xlabel('Speed [km/s]')
    hist_ax.set_ylabel('Frequency')


plot_hist(df['speed'], 20, vmin, vmax, cmap, hist_ax)
plt.tight_layout()
plt.savefig('fluxon-data/fss_theta_speed.png')
plt.show()