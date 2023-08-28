"""Plots the magnetic field strength and fluxon area
of the fluxons at the lower and upper boundaries

Parameters
----------
--cr: Carrington Rotation
--dat_dir: data directory
--show: (bool) show the plot
--batch: batch name
--file: file name
--nwant: number of fluxons to plot
"""

import os.path
import argparse
import matplotlib.pyplot as plt
import numpy as np

import py_plot_helper
from py_pipe_helper import (get_fixed_coords, load_fits_magnetogram, load_magnetogram_params)


print("\n\tPlotting Bmag...", end="")
# create the argument parser
parser = argparse.ArgumentParser(description=
                                'This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=0, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default=
                    '/Users/cgilbert/vscode/fluxons/fluxon-data', help='data directory')
parser.add_argument('--show', type=int, default=0)
parser.add_argument('--batch', type=str, default='default_batch')
parser.add_argument('--file', type=str, default="/Users/cgilbert/vscode/fluxons/fluxon-data/batches/default_batch/cr2193/wind/cr2193_f500_radial_bmag.dat")
parser.add_argument('--nwant', type=int, default=None, help='magnetogram file')

args = parser.parse_args()
batch = args.batch
(hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)
CR = args.cr or cr or 2183

filename = args.file or f'{args.dat_dir}/batches/{batch}/cr{CR}/wind/ \
                            cr{args.cr}_f{args.nwant}_radial_bmag.dat'

# Load the dat file
arr = np.loadtxt(filename).T
fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1 = arr
nfluxon = arr.shape[1]

# Convert coords to correct coords
ph0, th0 = get_fixed_coords(phi0, theta0)
ph1, th1 = get_fixed_coords(phi1, theta1)

# Do some data manipulation
br0_max = np.nanmax(br0) or 0.25
br1_max = np.nanmax(br1) or 0.25
ar0_max = np.nanmax(ar0) or 0.25
ar1_max = np.nanmax(ar1) or 0.25

skew = 5**2
power = 1
b0 = skew*(6*br0/br0_max)**power
b1 = skew*(4*br1/br1_max)**power
a0 = skew*(6*ar0/ar0_max)**power
a1 = skew*(4*ar1/ar1_max)**power

# Plot the Data
fig, (ax0, ax1) = plt.subplots(2)

magnet = load_fits_magnetogram()
ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
           extent=(0,2*np.pi,-1,1), aspect='auto')
ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower",
           extent=(0,2*np.pi,-1,1), aspect='auto')

sc00 = ax0.scatter(ph0, th0, c=br0, s = b0, cmap="winter",
                   alpha=0.75, label=r'B(1.0R$_\{odot}$)')
sc01 = ax0.scatter(ph1, th1, c=br1, s = b1, cmap="autumn",
                   alpha=0.75, label=r'B(21.5R$_\{odot}$)', marker='s')
sc10 = ax1.scatter(ph0, th0, c=ar0, s = a0, cmap="winter",
                   alpha=0.75, label=r'A(1.0R$_\{odot}$)')
sc11 = ax1.scatter(ph1, th1, c=ar1, s = a1, cmap="autumn",
                   alpha=0.75, label=r'A(21.5R$_\{odot}$)', marker='s')

cbar01 = fig.colorbar(sc01, ax=ax0)
cbar00 = fig.colorbar(sc00, ax=ax0)
cbar11 = fig.colorbar(sc11, ax=ax1)
cbar10 = fig.colorbar(sc10, ax=ax1)

cbar00.set_label(r"B(1.0R$_\{odot}$)  [Gauss] ")
cbar01.set_label(r"B(21.5R$_\{odot}$) [Gauss]")
cbar10.set_label(r"A(1.0R$_\{odot}$)  [m$^2$]")
cbar11.set_label(r"A(21.5R$_\{odot}$) [m$^2$]")

ax0.set_title(F"Fluxon Magnetic Field Strength for CR {args.cr}")
ax1.set_title(F"Fluxon Area for CR {args.cr}")

for ax in (ax0, ax1):
    ax.set_xlabel('Longitude (Radians)')
    ax.set_ylabel('Sine latitude')
    ax.set_ylim((-1.5,1.1))
    ax.axhline(-1, c='lightgrey', zorder=-10)
    ax.axhline( 1, c='lightgrey', zorder=-10)
    ax.axvline(0, c='lightgrey', zorder=-10)
    ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
    ax.legend(loc="lower right")

fig.set_size_inches((8,8))
plt.tight_layout()

imagename = os.path.basename(filename.replace(".dat", ".png"))
imagedir = os.path.dirname(os.path.dirname(os.path.dirname(filename)))
bdir = os.path.join(imagedir, "imgs", "bmag")
if not os.path.exists(bdir):
    os.makedirs(bdir)
pngname = os.path.join(bdir, imagename)

plt.savefig(pngname)
if args.show:
    plt.show()
plt.close(fig)
print("Done!\n")
print("Saved to ", pngname, "\n")
