# read_fr.py
# Script to parse fluxon output file of the field expansion factor


# file output format: fps, phb, thb, phe, the, vrv, vre, frb, fre, fre2
# fluxon id, beginning coords, end coords, beginning velocity, end velocity, beginning expansion, end expansion

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
from magnetoget import load_fits_magnetogram

default_cr = 2163

# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=default_cr, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/Fluxon-Scripts-Gilly', help='data directory')
parser.add_argument('--show', type=int, default=1)
args = parser.parse_args()
filename = f'{args.dat_dir}/fluxon/cr{args.cr}/wind/radial_wind.dat'

# Load the dat file
arr = np.loadtxt(filename).T
try:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
except ValueError:
    fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
nfluxon = arr.shape[1]

# Convert coords to correct coords
ph0, th0 = phi0+np.pi, -theta0+(np.pi/2)
ph1, th1 = phi1+np.pi, -theta1+(np.pi/2)


# Do some data manipulation

#remove strongest
np.sort(vel0.flatten())


vel0[vel0<0.]= np.mean(vel0[vel0>0.])
vel1[vel1<0.]= np.mean(vel1[vel1>0.])
fr0[fr0<0.]  = np.nan
fr1[fr1<0.]  = np.nan

arg0 = np.argwhere(vel0[vel0<0.])
arg1 = np.argwhere(vel1[vel1<0.])
arg2 = np.argwhere(fr0[fr0<0.]  )
arg3 = np.argwhere(fr1[fr1<0.]  )

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
fig, (ax0, ax1, ax2) = plt.subplots(3, sharex = True, sharey= False)

magnet = load_fits_magnetogram()
ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')

sc00 = ax0.scatter(ph0, th0, c= vel0, s=v0, alpha=0.75, label='V(1.0R$_\odot$)'  , cmap="winter",  )
sc01 = ax0.scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
# sc01 = ax1.scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
# sc10 = ax1.scatter(ph0, th0, c= fr0 , s=f0, alpha=0.75, label='Fr(1.0R$_\odot$)' , cmap="winter", )
# sc11 = ax1.scatter(ph1, th1, c= fr1 , s=f1, alpha=0.75, label='Fr(21.5R$_\odot$)', cmap="autumn", marker='s')


from scipy.interpolate import griddata

# ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto', zorder=-10)
sc00 = ax1.scatter(ph0, th0, c= vel0, s=v0, alpha=1, label='V(1.0R$_\odot$)'  , cmap="winter",  edgecolors='k')
sc01 = ax2.scatter(ph1, th1, c= vel1, s=v1, alpha=1, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s',  edgecolors='k')


# Create a grid for interpolation
Ny, Nx = magnet.shape
grid_x, grid_y = np.linspace(0, 2*np.pi, Nx), np.linspace(-1, 1,  Ny)
grid_x, grid_y = np.meshgrid(grid_x, grid_y)

# Interpolate values on the grid
points0 = [(ph, th) for ph, th in zip(ph0,th0)]
points1 = [(ph, th) for ph, th in zip(ph1,th1)]
grid_z0 = griddata(points0, vel0, (grid_x, grid_y), method='nearest', fill_value=0)
grid_z1 = griddata(points1, vel1, (grid_x, grid_y), method='nearest', fill_value=0)

# Create a contour plot
# fig, ax = plt.subplots()
contour0 = ax1.contourf(grid_x, grid_y, grid_z0, zorder=0, alpha = 1, cmap="winter")
contour1 = ax2.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha = 1, cmap="autumn")

# Add a colorbar
cbar0 = fig.colorbar(contour0, ax=ax1)
cbar1 = fig.colorbar(contour1, ax=ax2)
cbar0.set_label('Interpolated Wind [km/s]')
cbar1.set_label('Interpolated Wind [km/s]')

# ax0.set_



# cbar01 = fig.colorbar(sc01, ax=ax0)
# cbar00 = fig.colorbar(sc00, ax=ax0)
# cbar11 = fig.colorbar(sc11, ax=ax1)
# cbar10 = fig.colorbar(sc10, ax=ax1)

# cbar00.set_label("V(1.0R$_\odot$)  [km/s]")
# cbar01.set_label("V(21.5R$_\odot$) [km/s]")
# cbar10.set_label("Fr(1.0R$_\odot$)  ")
# cbar11.set_label("Fr(21.5R$_\odot$) ")


fig.suptitle(F"Fluxon Wind Strength for CR {args.cr}")
# ax1.set_title(F"Fluxon Expansion for CR {args.cr}")
for ax in [ax0, ax1, ax2]:
    ax.set_xlabel('Longitude (Radians)')
    ax.set_ylabel('Sine latitude')
    ax.set_ylim((-1.5,1.1))
    ax.axhline(-1, c='lightgrey', zorder=-10)
    ax.axhline( 1, c='lightgrey', zorder=-10)
    ax.axvline(0, c='lightgrey', zorder=-10)
    ax.axvline(2*np.pi, c='lightgrey', zorder=-10)

ax0.legend(loc="lower right")
ax1.set_ylim((-1, 1))
ax2.set_ylim((-1, 1))

fig.set_size_inches((8,8))
plt.tight_layout()
pngname = filename.replace(".dat", "_map.png")
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