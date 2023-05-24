# read_fr.py
# Script to parse fluxon output file of the field expansion factor


# file output format: fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1
# fluxon id, beginning coords, end coords, beginning mag, end mag, area maybe?

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse

default_cr = 2160

# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=default_cr, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/Fluxon-Scripts-Gilly', help='data directory')
parser.add_argument('--show', type=int, default=0)
args = parser.parse_args()
filename = f'{args.dat_dir}/fluxon/cr{args.cr}/wind/radial_bmag.dat'

# Load the dat file
arr = np.loadtxt(filename).T
fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1 = arr
nfluxon = arr.shape[1]

fig, (ax0, ax1) = plt.subplots(2)


br0_max = np.nanmax(br0) or 0.25
br1_max = np.nanmax(br1) or 0.25
ar0_max = np.nanmax(ar0) or 0.25
ar1_max = np.nanmax(ar1) or 0.25


ax0.set_title(F"Fluxon Magnetic Field Strength for CR {args.cr}")
# ax0.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*br0/br0_max)**3, alpha=0.75, label='B(1.0R$_\odot$)')
# ax0.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*br1/br1_max)**3, alpha=0.75, label='B(21.5R$_\odot$)')
# ax1.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*ar0/ar0_max)**3, alpha=0.75, label='A(1.0R$_\odot$)', marker='s')
# ax1.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*ar1/ar1_max)**3, alpha=0.75, label='A(21.5R$_\odot$)', marker='s')

skew = 2
b0 = (skew+5*br0/br0_max)**3
b1 = (skew+5*br1/br1_max)**3
a0 = (skew+5*ar0/ar0_max)**3
a1 = (skew+5*ar1/ar1_max)**3

ax0.scatter(phi0+np.pi, -(theta0-(np.pi/2)), s=b0, alpha=0.75, label='B(1.0R$_\odot$)')
ax0.scatter(phi1+np.pi, -(theta1-(np.pi/2)), s=b1, alpha=0.75, label='B(21.5R$_\odot$)')
ax1.scatter(phi0+np.pi, -(theta0-(np.pi/2)), s=a0, alpha=0.75, label='A(1.0R$_\odot$)', marker='s')
ax1.scatter(phi1+np.pi, -(theta1-(np.pi/2)), s=a1, alpha=0.75, label='A(21.5R$_\odot$)', marker='s')


ax1.set_title(F"Fluxon (Area?) for CR {args.cr}")
for ax in (ax0, ax1):
    ax.set_xlabel('Longitude (Radians)')
    ax.set_ylabel('Sine latitude')
    ax.set_ylim((-1.1,1.1))
    ax.axhline(-1, c='lightgrey', zorder=-10)
    ax.axhline( 1, c='lightgrey', zorder=-10)
    ax.axvline(0, c='lightgrey', zorder=-10)
    ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
    ax.legend(loc="lower right")

fig.set_size_inches((6,8))
plt.tight_layout()
pngname = filename.replace(".dat", ".png")
plt.savefig(pngname)
if args.show:
    plt.show()
plt.close(fig)