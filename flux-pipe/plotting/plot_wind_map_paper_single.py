"""_summary_

Returns
-------
_type_
    _description_

Raises
------
to
    _description_
"""
import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
from py_pipe_helper import (load_fits_magnetogram, load_magnetogram_params, get_fixed_coords)
import os.path as path
from plot_wind_map_paper_brief import remove_outliers
from scipy.interpolate import griddata
from py_plot_helper import scale_data


# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=2160, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxons/fluxon-data', help='data directory')
parser.add_argument('--show', type=int, default=1)
parser.add_argument('--interp', type=str, default="linear")
parser.add_argument('--nact', type=int, default=0)
parser.add_argument('--batch', type=str, default="fluxon_paperfigs", help='select the batch name')
args = parser.parse_args()
batch = args.batch
interp = args.interp

# Load the magnetogram parameters
(hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)
CR = args.cr or cr or 2183
print(f"plotting CR{CR}...", end="\n" if __name__=="__main__" else "")

# Load the wind file
dat_file = f'{args.dat_dir}/{batch}/cr{CR}/wind/radial_wind_f{args.nact}.dat'
arr = np.loadtxt(dat_file).T
try:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
except ValueError:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
nfluxon = arr.shape[1]

# Set the output file names
main_file =  f'{args.dat_dir}/{batch}/cr{CR}/wind/radial_wind_map_{len(vel1)}_{interp}_paper.png'
outer_file = f"{args.dat_dir}/{batch}/imgs/windmap_2/cr{CR}_radial_wind_map_ou{len(vel1)}_paper.png"
hist_file =  f"{args.dat_dir}/{batch}/cr{CR}/wind/vel_hist_{len(vel1)}_paper.png"

import os
if not path.exists(os.path.dirname(outer_file)):
    os.makedirs(os.path.dirname(outer_file))

# Convert coords to correct coords
ph0, th0 = get_fixed_coords(phi0, theta0)
ph1, th1 = get_fixed_coords(phi1, theta1)

vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 2, 1)
vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)



v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

fig, ax = plt.subplots(3)

do_hist = False
if do_hist:
## The Histogram Plot
    mean1 = np.mean(vel1_clean)
    std1  =  np.std(vel1_clean)
    ax[3].hist(vel1_clean[vel1_clean>=0], bins=10, color='sandybrown')
    ax[3].axvline(mean1, color='k', linestyle='dashed', linewidth=1,
                  label="mean: {mean1:.2f}")
    ax[3].axvline(mean1+std1, color='k', linestyle=':', linewidth=1,
                  label="std: {std1:.2f}")
    ax[3].axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
    ax[3].legend()
    ax[3].set_xlabel("Velocity (km/s)")
    ax[3].set_ylabel("Number of Fluxons")


## The Magnetogram Plot
magnet = load_fits_magnetogram(batch=batch, bo=3, bn=2)
# find the max and min of the magnetogram plot for use in setting the colormap,
sigma = 3
mmean = np.nanmean(magnet)
msig = np.nanstd(magnet)
mvmin = mmean - sigma*msig
mvmax = mmean + sigma*msig

# Plot the magnetogram
magimg = ax[0].imshow(magnet, cmap='gray', interpolation=None,
                      origin="lower", extent=(0,2*np.pi,-1,1),
                    aspect='auto', vmin=mvmin, vmax=mvmax)
# Scatter the fluxon data
ph000 = (ph0_clean+np.pi)%(2*np.pi)
ph111 = (ph1_clean+np.pi)%(2*np.pi)
ph000b = (ph0b+np.pi)%(2*np.pi)
ph111b = (ph1b+np.pi)%(2*np.pi)

sc00 = ax[0].scatter(ph000, th0_clean, c= vel0_clean, s=v0, alpha=0.5,
                     label=r'V(1.0R$_\odot$)'  , cmap="winter",  )
sc01 = ax[0].scatter(ph111, th1_clean, c= vel1_clean, s=v1, alpha=0.5,
                     label=r'V(21.5R$_\odot$)' , cmap="autumn",  marker='s')

#Scatter the Outliers
ax[0].scatter(ph000b, th0b, c= v0bs, alpha=0.6,
              label=r'V(1.0R$_\odot$)', cmap="winter", marker='X',  s=50, edgecolors='k')
ax[0].scatter(ph111b, th1b, c= v1bs, alpha=0.6,
              label=r'V(1.5R$_\odot$)', cmap="autumn", marker='X',  s=50, edgecolors='k')
# ax[1].scatter(ph000b, th0b, c= outlier_V0_scaled, alpha=0.95, label='V(1.0R$_\odot$)',
# cmap="winter", marker='X', s=50, edgecolors='k')
ax[1].scatter(ph111b, th1b, c= v1bs, alpha=0.95,
              label=r'V(1.5R$_\odot$)', cmap="autumn", marker='X', s=50, edgecolors='k')


## Plot the Interp. data
# Define the x-boundaries of the domain
x_min = 0
x_max = 2*np.pi

# Wrap the data around the edges of the domain
ph0 = ph000
th0 = th0_clean
ph1 = ph111
th1 = th1_clean
vel0 = vel0_clean
vel1 = vel1_clean


ph0_wrapped = np.concatenate((ph0 - x_max, ph0, ph0 + x_max))
th0_wrapped = np.concatenate((th0, th0, th0))
vel0_wrapped = np.concatenate((vel0, vel0, vel0))

ph1_wrapped = np.concatenate((ph1 - x_max, ph1, ph1 + x_max))
th1_wrapped = np.concatenate((th1, th1, th1))
vel1_wrapped = np.concatenate((vel1, vel1, vel1))

# Create a grid for interpolation
Ny, Nx = magnet.shape
grid_x, grid_y = np.linspace(x_min, x_max, Nx, endpoint=False), np.linspace(-1, 1, Ny)
grid_x, grid_y = np.meshgrid(grid_x, grid_y)

# Interpolate values on the grid
points0 = [(ph, th) for ph, th in zip(ph0_wrapped, th0_wrapped)]
points1 = [(ph, th) for ph, th in zip(ph1_wrapped, th1_wrapped)]
grid_z0 = griddata(points0, vel0_wrapped, (grid_x, grid_y), method=interp, fill_value=0)
grid_z1 = griddata(points1, vel1_wrapped, (grid_x, grid_y), method=interp, fill_value=0)

# Create a scatter plot
sc01 = ax[1].scatter(ph1, th1, c=vel1, s=v1, alpha=1,
                     label=r'V(21.5R$_\odot$)', cmap="autumn", marker='s', edgecolors='k')

# Create a contour plot
# contour0 = ax[1].contourf(grid_x, grid_y, grid_z0, zorder=0, alpha=1, cmap="winter")
contour1 = ax[1].contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap="autumn")



# Add a colorbar
# cbar0 = fig.colorbar(contour0, ax=ax[1])
cbar1 = fig.colorbar(contour1, ax=ax[1], pad=0.03)
# cbar0.set_label('Interp. Wind [km/s]\nV( 1.0R$_\odot$)')
cbar1.set_label('Interp. Wind Speed \n[km/s]') #\nV(21.5R$_\odot$)


# ax[0].set_


cbar01 = fig.colorbar(magimg, ax=ax[0], extend='both', pad=0.03)
cbar01.set_label("Mag Field Anomaly [Gauss]")

nx = 20
sy, sx = grid_x.shape
ratio = float(sy)/float(sx)
gridsize = (int(nx*1.2), int(np.round(nx*ratio)))
hex1 = ax[2].hexbin(grid_x.flatten(), grid_y.flatten(), C=grid_z1.flatten(),
        gridsize=gridsize, cmap='autumn', vmin=np.min(grid_z1), vmax=np.max(grid_z1))

cbar2 = fig.colorbar(hex1, ax=ax[2], pad=0.03)
cbar2.set_label("Interp. Wind Speed \n[km/s]")

fig.suptitle(F"Fluxon Wind Strength for CR {CR}: Open Fluxons={len(v1)}")
# ax[1].set_title(F"Fluxon Expansion for CR {args.cr}")
for this_ax in [ax[0], ax[1], ax[2]]:
    this_ax.set_ylabel('Sine latitude')
    this_ax.set_ylim((-1.0,1.0))
    this_ax.set_aspect('equal')
    this_ax.axhline(-1, c='lightgrey', zorder=-10)
    this_ax.axhline( 1, c='lightgrey', zorder=-10)
    this_ax.axvline(0, c='lightgrey', zorder=-10)
    this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)

ax[2].set_xlabel('Longitude (Radians)')
ax[1].set_xlim((0, 2*np.pi))
ax[2].set_xlim((0, 2*np.pi))
# ax[0].legend(loc="lower right")
ax[0].set_ylim((-1, 1))
ax[1].set_ylim((-1, 1))


fig.set_size_inches((8.5,8))
# plt.show()

plt.tight_layout()

# Save the Figures
print("Saving figures...")
print(main_file)
# plt.savefig(main_file, dpi=200)
plt.savefig(main_file.replace(".png", ".pdf"), dpi=200)

print(outer_file)
# plt.savefig(outer_file, dpi=200)
plt.savefig(outer_file.replace(".png", ".pdf"), dpi=200)

plt.close(fig)
