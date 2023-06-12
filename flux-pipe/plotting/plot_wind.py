# read_fr.py
# Script to parse fluxon output file of the field expansion factor


# file output format: fps, phb, thb, phe, the, vrv, vre, frb, fre, fre2
# fluxon id, beginning coords, end coords, beginning velocity, end velocity, beginning expansion, end expansion

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
from py_pipe_helper import load_fits_magnetogram

default_cr = 2163

# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=default_cr, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
parser.add_argument('--batch', type=str, default="fluxon", help='select the batch name')
parser.add_argument('--show', type=int, default=1)
args = parser.parse_args()
batch = args.batch
filename = f'{args.dat_dir}/{batch}/cr{args.cr}/wind/radial_wind.dat'

# Load the dat file
arr = np.loadtxt(filename).T
try:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
except ValueError:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
nfluxon = arr.shape[1]


# Convert coords to correct coords
from py_pipe_helper import get_fixed_coords
ph0, th0 = get_fixed_coords(phi0, theta0)
ph1, th1 = get_fixed_coords(phi1, theta1)

# ph0, th0 = phi0+np.pi, np.sin(-(theta0-(np.pi/2)))
# ph1, th1 = phi1+np.pi, np.sin(-(theta1-(np.pi/2)))


# Do some data manipulation

vel0[vel0<0.]=np.nan
vel1[vel1<0.]=np.nan
fr0[fr0<0.]  =np.nan
fr1[fr1<0.]  =np.nan

vel0_max = np.nanmax(np.abs(vel0)) 
vel1_max = np.nanmax(np.abs(vel1)) 
fr0_max  = np.nanmax(np.abs(fr0))  
fr1_max  = np.nanmax(np.abs(fr1))  

vel0_min = np.nanmin(np.abs(vel0)) 
vel1_min = np.nanmin(np.abs(vel1)) 
fr0_min  = np.nanmin(np.abs(fr0))  
fr1_min  = np.nanmin(np.abs(fr1)) 

scale = 15**2
power = 1

v0 = scale * ((np.abs(vel0) - vel0_min) / (vel0_max-vel0_min))**power
v1 = scale * ((np.abs(vel1) - vel1_min) / (vel1_max-vel1_min))**power
f0 = scale * ((np.abs(fr0 ) - fr0_min ) / (fr0_max -fr0_min ))**power
f1 = scale * ((np.abs(fr1 ) - fr1_min ) / (fr1_max -fr1_min ))**power

# Plot the Data
fig, (ax0, ax1) = plt.subplots(2)

magnet = load_fits_magnetogram()
ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')
ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')


sc00 = ax0.scatter(ph0, th0, c= vel0, s=v0, alpha=0.75, label='V(1.0R$_\odot$)'  , cmap="winter",  )
sc01 = ax0.scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
sc10 = ax1.scatter(ph0, th0, c= fr0 , s=f0, alpha=0.75, label='Fr(1.0R$_\odot$)' , cmap="winter", )
sc11 = ax1.scatter(ph1, th1, c= fr1 , s=f1, alpha=0.75, label='Fr(21.5R$_\odot$)', cmap="autumn", marker='s')

cbar01 = fig.colorbar(sc01, ax=ax0)
cbar00 = fig.colorbar(sc00, ax=ax0)
cbar11 = fig.colorbar(sc11, ax=ax1)
cbar10 = fig.colorbar(sc10, ax=ax1)

cbar00.set_label("V(1.0R$_\odot$)  [km/s]")
cbar01.set_label("V(21.5R$_\odot$) [km/s]")
cbar10.set_label("Fr(1.0R$_\odot$)  ")
cbar11.set_label("Fr(21.5R$_\odot$) ")


ax0.set_title(F"Fluxon Wind Strength for CR {args.cr}")
ax1.set_title(F"Fluxon Expansion for CR {args.cr}")
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
pngname = filename.replace(".dat", ".png")
plt.savefig(pngname)
# if args.show:
#     plt.show()
plt.close(fig)




# ax0.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*vel0/vel0.max())**3, alpha=0.75, label='V(1.0R$_\odot$)')
# ax0.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*vel1/vel1.max())**3, alpha=0.75, label='V(21.5R$_\odot$)')
# ax1.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*fr0/fr0.max())**3, alpha=0.75, label='Fr(1.0R$_\odot$)', marker='s')
# ax1.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*fr1/fr1.max())**3, alpha=0.75, label='Fr(21.5R$_\odot$)', marker='s')


 

# v0 = (1+5*np.abs(vel0)/vel0_max)**3
# v1 = (1+5*np.abs(vel1)/vel1_max)**3
# f0 = (1+5*np.abs(fr0 )/fr0_max )**3
# f1 = (1+5*np.abs(fr1 )/fr1_max )**3