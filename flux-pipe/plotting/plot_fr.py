# read_fr.py
# Script to parse fluxon output file of the field expansion factor

# file output format: fid, x, y, z, r, th, ph, a

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse

default_cr = 2163

# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=default_cr, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
parser.add_argument('--batch', type=str, default="fluxon", help='select the batch name')
parser.add_argument('--show', type=int, default=0)
args = parser.parse_args()
batch = args.batch
filename = f'{args.dat_dir}/batches/{batch}/cr{args.cr}/wind/radial_fr.dat'
# print(filename)

# Load the dat file
arr = np.loadtxt(filename).T
nfluxon = int(arr[0,:].max())

# Initialize Empty Arrays
th0 = np.zeros(nfluxon)
ph0 = np.zeros(nfluxon)
a0  = np.zeros(nfluxon)
th1 = np.zeros(nfluxon)
ph1 = np.zeros(nfluxon)
a1  = np.zeros(nfluxon)

# Fill the arrays
for i in np.arange(0, nfluxon):
    floc = np.where(arr[0,:] == i)[0]

    th0[i] = -1 * arr[5, floc[0]] + np.pi/2
    th1[i] = -1 * arr[5, floc[-1]] + np.pi/2

    ph0[i] = arr[6, floc[0]] + np.pi
    ph1[i] = arr[6, floc[-1]] + np.pi

    a0[i] = arr[7, floc[0]]
    a1[i] = arr[7, floc[-1]]

# a1[np.where(a1 > 1e21)] = 0


## Plot things!

fig, (ax1, ax2) = plt.subplots(2)
# TODO make sure the X axis is correct
shift = -1.0
ax1.set_title(F"Fluxon Expansion Factors for CR {args.cr}")
RS = 696340000 #meters
for i in np.arange(0, nfluxon):
    floc = np.where(arr[0,:] == i)[0]
    ax1.plot((shift + arr[4,floc]/RS), arr[7,floc])
r1 = 1.
r2 = 21.5
ax1.axvline(r1, ls=":", c='lightgrey')
ax1.axvline(r2, ls=":", c='lightgrey')
ax1.set_xlabel('Heliocentric Radius [R$_\odot$]')
ax1.set_ylabel('Cross Sectional Area [m$^2$]')
ax1.set_yscale('log')


ax2.set_title("Fluxon Area at Low (1R$_\odot$) and High (21.5R$_\odot$) Radii")

a0_max = np.nanmax(np.abs(a0)) or 0.25
a1_max = np.nanmax(np.abs(a1)) or 0.25
skew = 2

aa0 = (skew+5*np.abs(a0)/a0_max)**3
aa1 = (skew+5*np.abs(a1)/a1_max)**3

ax2.scatter(ph0, np.sin(th0), s = aa0 , label='A(1.0R$_\odot$),', alpha=0.75, marker='s')
ax2.scatter(ph1, np.sin(th1), s = aa1 , label='A(21.5R$_\odot$)', alpha=0.75, marker='s')
ax2.set_xlabel('Longitude (Radians)')
ax2.set_ylabel('Sine latitude')
ax2.set_ylim((-1.1,1.1))
ax2.axhline(-1, c='lightgrey', zorder=-10)
ax2.axhline( 1, c='lightgrey', zorder=-10)
ax2.axvline(0, c='lightgrey', zorder=-10)
ax2.axvline(2*np.pi, c='lightgrey', zorder=-10)
ax2.legend(loc="upper right")

fig.set_size_inches((6,8))
plt.tight_layout()
pngname = filename.replace(".dat", ".png")
plt.savefig(pngname)
if args.show:
    plt.show()
plt.close(fig)