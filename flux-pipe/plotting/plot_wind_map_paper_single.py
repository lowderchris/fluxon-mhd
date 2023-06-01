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
import os.path as path
from magnetoget import load_magnetogram_params
from plot_wind_map_paper_brief import remove_outliers


# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=2160, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/Fluxon-Scripts-Gilly', help='data directory')
parser.add_argument('--show', type=int, default=1)
parser.add_argument('--interp', type=str, default="linear")
parser.add_argument('--nact', type=int, default=0)
parser.add_argument('--batch', type=str, default="fluxon_paperfigs", help='select the batch name')
args = parser.parse_args()
batch = args.batch
interp = args.interp
# batch = "fluxon_HQ_2000_2"


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
ph0, th0 = phi0+np.pi, -theta0+(np.pi/2)
ph1, th1 = phi1+np.pi, -theta1+(np.pi/2)


# # Remove outliers from the dataset recursively
# def remove_outliers(data, ph, th, threshold=3, n_times= 1):

#     # Get the Good Points
#     mean = np.mean(data[data>0])
#     std = np.std(data[data>0])
#     data_good, inds_good = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data) if (mean - threshold * std < x < mean + threshold * std) and x>0]))
#     ph_good, th_good = ph[inds_good], th[inds_good]
#     data_good = np.asarray(data_good)

#     # Get the Bad Points
#     inds_bad = [i for i in range(len(data)) if i not in inds_good]
#     ph_bad, th_bad = ph[inds_bad], th[inds_bad]
#     data_bad = data[inds_bad]

#     # Return or Recurse
#     if n_times <= 1:
#         return data_good, ph_good, th_good, data_bad, ph_bad, th_bad
#     else:
#         data_good_2, ph_good_2, th_good_2, data_bad_2, ph_bad_2, th_bad_2 = remove_outliers(data_good, ph_good, th_good, threshold, n_times-1)
#         bad_data_all = np.concatenate((data_bad, data_bad_2))
#         bad_ph_all   = np.concatenate((  ph_bad,   ph_bad_2))
#         bad_th_all   = np.concatenate((  th_bad,   th_bad_2))
#         return data_good_2, ph_good_2, th_good_2, bad_data_all, bad_ph_all, bad_th_all

# vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 1)
# vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2)

vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 2, 1)
vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)


def scale_data(vel0_clean, vel1_clean, outlier_V0, outlier_V1, scale=15**2, power=1):
    """Scale the data between 0 and 1, then raise to a power, then scale by a factor""" 
    vel0_max = np.nanmax(vel0_clean)
    vel1_max = np.nanmax(vel1_clean)
    vel0_min = np.nanmin(vel0_clean)
    vel1_min = np.nanmin(vel1_clean)

    v0 = scale * ((np.abs(vel0_clean) - vel0_min) / (vel0_max-vel0_min))**power
    v1 = scale * ((np.abs(vel1_clean) - vel1_min) / (vel1_max-vel1_min))**power

    outlier_V0_scaled = scale * ((np.abs(outlier_V0) - vel0_min) / (vel0_max-vel0_min))**power
    outlier_V1_scaled = scale * ((np.abs(outlier_V1) - vel1_min) / (vel1_max-vel1_min))**power
    
    return v0, v1, outlier_V0_scaled, outlier_V1_scaled

v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)




import matplotlib.gridspec as gridspec

fig, ax = plt.subplots(3)
## The Histogram Plot

# mean1 = np.mean(vel1_clean)
# std1  =  np.std(vel1_clean)
# ax[3].hist(vel1_clean[vel1_clean>=0], bins=10, color='sandybrown')
# ax[3].axvline(mean1, color='k', linestyle='dashed', linewidth=1, label="mean: {:.2f}".format(mean1))
# ax[3].axvline(mean1+std1, color='k', linestyle=':', linewidth=1, label="std: {:.2f}".format(std1))
# ax[3].axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
# ax[3].legend()
# # ax[3].set_title(f"Solar Wind Velocity Distribution: CR{CR}, Open Fluxons: {len(vel1)}", fontsize=16)
# ax[3].set_xlabel("Velocity (km/s)")
# ax[3].set_ylabel(F"Number of Fluxons")


## The Magnetogram Plot
magnet = load_fits_magnetogram(batch=batch, bo=3, bn=2)
# find the max and min of the magnetogram plot for use in setting the colormap, 
sigma = 3
mmean = np.nanmean(magnet)
msig = np.nanstd(magnet)
mvmin = mmean - sigma*msig
mvmax = mmean + sigma*msig

# Plot the magnetogram
magimg = ax[0].imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), 
                    aspect='auto', vmin=mvmin, vmax=mvmax)
# Scatter the fluxon data
ph000 = (ph0_clean+np.pi)%(2*np.pi)
ph111 = (ph1_clean+np.pi)%(2*np.pi)
ph000b = (ph0b+np.pi)%(2*np.pi)
ph111b = (ph1b+np.pi)%(2*np.pi)

sc00 = ax[0].scatter(ph000, th0_clean, c= vel0_clean, s=v0, alpha=0.5, label='V(1.0R$_\odot$)'  , cmap="winter",  )
sc01 = ax[0].scatter(ph111, th1_clean, c= vel1_clean, s=v1, alpha=0.5, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
# sc01 = ax[1].scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
# sc10 = ax[1].scatter(ph0, th0, c= fr0 , s=f0, alpha=0.75, label='Fr(1.0R$_\odot$)' , cmap="winter", )
# sc11 = ax[1].scatter(ph1, th1, c= fr1 , s=f1, alpha=0.75, label='Fr(21.5R$_\odot$)', cmap="autumn", marker='s')

#Scatter the Outliers
ax[0].scatter(ph000b, th0b, c= v0bs, alpha=0.6, label='V(1.0R$_\odot$)', cmap="winter", marker='X',  s=50, edgecolors='k')
ax[0].scatter(ph111b, th1b, c= v1bs, alpha=0.6, label='V(1.5R$_\odot$)', cmap="autumn", marker='X',  s=50, edgecolors='k')
# ax[1].scatter(ph000b, th0b, c= outlier_V0_scaled, alpha=0.95, label='V(1.0R$_\odot$)', cmap="winter", marker='X', s=50, edgecolors='k')
ax[1].scatter(ph111b, th1b, c= v1bs, alpha=0.95, label='V(1.5R$_\odot$)', cmap="autumn", marker='X', s=50, edgecolors='k')



from scipy.interpolate import griddata

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
# sc00 = ax[1].scatter(ph0, th0, c=vel0, s=v0, alpha=1, label='V(1.0R$_\odot$)', cmap="winter", edgecolors='k'             )
sc01 = ax[1].scatter(ph1, th1, c=vel1, s=v1, alpha=1, label='V(21.5R$_\odot$)', cmap="autumn", marker='s', edgecolors='k')

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

# cbar01 = fig.colorbar(sc01, ax=ax[0])
# cbar00 = fig.colorbar(sc00, ax=ax[0])
# cbar11 = fig.colorbar(sc11, ax=ax[1])
# cbar10 = fig.colorbar(sc10, ax=ax[1])

# cbar00.set_label("V(1.0R$_\odot$)  [km/s]")
# cbar01.set_label("V(21.5R$_\odot$) [km/s]")
# cbar10.set_label("Fr(1.0R$_\odot$)  ")
# cbar11.set_label("Fr(21.5R$_\odot$) ")



# axarr[2].set_title("hexbin")

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
# ax[3].set_ylim((-1, 1))


        

# ax[3].set_position(ax[2].get_position())

fig.set_size_inches((8.5,8))
# plt.show()

plt.tight_layout()

# plt.subplots_adjust(
# top=0.955,
# bottom=0.046,
# left=0.095,
# right=0.98,
# hspace=0.242,
# wspace=0.15
# )


# Save the Figures
print("Saving figures...")
print(main_file)
# plt.savefig(main_file, dpi=200)
plt.savefig(main_file.replace(".png", ".pdf"), dpi=200)

print(outer_file)
# plt.savefig(outer_file, dpi=200)
plt.savefig(outer_file.replace(".png", ".pdf"), dpi=200)

plt.close(fig)


# def plot_three(grid_x, grid_y, grid_z, save_path):
#     fig, axarr = plt.subplots(3, sharex="all", sharey="all", figsize=(10,10))
#     # the imshow extent will let it plot on the same space as the other two plots, using grid_x and grid_y
#     extent = (np.min(grid_x), np.max(grid_x), np.min(grid_y), np.max(grid_y))
#     axarr[0].set_title("imshow")
#     axarr[0].imshow(grid_z, cmap ="autumn", extent=extent, interpolation='None', origin='lower')
#     axarr[1].set_title("contourf")
#     axarr[1].contourf(grid_x, grid_y, grid_z, zorder=0, alpha=1, cmap="autumn")
#     axarr[2].set_title("hexbin")

#     nx = 12
#     sy, sx = grid_x.shape
#     ratio = float(sy)/float(sx)
#     gridsize = (nx, int(np.round(nx*ratio)))

#     axarr[2].hexbin(grid_x.flatten(), grid_y.flatten(), C=grid_z.flatten(), gridsize=gridsize, cmap='autumn')
    

#     plt.tight_layout()
#     plt.savefig(save_path)
#     plt.close(fig)
#     # plt.show()


# ax[0].scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*vel0/vel0.max())**3, alpha=0.75, label='V(1.0R$_\odot$)')
# ax[0].scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*vel1/vel1.max())**3, alpha=0.75, label='V(21.5R$_\odot$)')
# ax[1].scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*fr0/fr0.max())**3, alpha=0.75, label='Fr(1.0R$_\odot$)', marker='s')
# ax[1].scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*fr1/fr1.max())**3, alpha=0.75, label='Fr(21.5R$_\odot$)', marker='s')

# f0 = scale * ((np.abs(fr0 ) - fr0_min ) / (fr0_max -fr0_min ))**power
# f1 = scale * ((np.abs(fr1 ) - fr1_min ) / (fr1_max -fr1_min ))**power

 

# v0 = (1+5*np.abs(vel0)/vel0_max)**3
# v1 = (1+5*np.abs(vel1)/vel1_max)**3
# f0 = (1+5*np.abs(fr0 )/fr0_max )**3
# f1 = (1+5*np.abs(fr1 )/fr1_max )**3



# # Find the positive mean
# v0mean = np.nanmean(vel0[vel0>0.])
# v1mean = np.nanmean(vel1[vel1>0.])
# v0std = np.nanstd(vel0[vel0>0.])
# v1std = np.nanstd(vel1[vel1>0.])

# # Set the negative values to the mean
# vel0[vel0<0.]= v0mean
# vel1[vel1<0.]= v1mean
# fr0_max  = np.nanmax(np.abs(fr0))  
# fr1_max  = np.nanmax(np.abs(fr1))  
# fr0_min  = np.nanmin(np.abs(fr0))  
# fr1_min  = np.nanmin(np.abs(fr1)) 

# # Set the large values to less large
# sig=2
# v0big = (v0mean+sig*v0std)
# v1big = (v1mean+sig*v1std)
# vel0[vel0>v0big] = v0big
# vel1[vel1>v1big] = v1big



#############################################
#  Manipulate Data - More Art than Science  #
#############################################


# # import numpy as np
# vel0_clean, mean0, std0, bad_inds0 = remove_outliers(vel0)
# vel1_clean, mean1, std1, bad_inds1 = remove_outliers(vel1)
# vel1_clean, mean1, std1, bad_inds1b = remove_outliers(vel1_clean)

# # # Convert vel0 and vel1 to lists
# vel0_list = list(vel0_clean)
# vel1_list = list(vel1_clean)
# ph0_list = list(ph0)
# th0_list = list(th0)
# ph1_list = list(ph1)
# th1_list = list(th1)

# # Determine the indices of outliers in vel0 and vel1
# vel0_mean = np.mean(vel0[vel0>0])
# vel1_mean = np.mean(vel1[vel1>0])
# vel0_std = np.std(vel0[vel0>0])
# vel1_std = np.std(vel1[vel1>0])
# threshold = 2.5

# # This is a list of indices of outliers
# vel0_indices = [idx for idx, val in enumerate(vel0_list) if np.abs(val - vel0_mean) / vel0_std > threshold or val<0]
# vel1_indices = [idx for idx, val in enumerate(vel1_list) if np.abs(val - vel1_mean) / vel1_std > threshold or val<0]

# # Remove outliers from vel0 and vel1

# outlier_V0 = [];
# outPh0 = [];
# outTh0 = [];
# outlier_V1 = [];
# outPh1 = [];
# outTh1 = [];



# for idx in sorted(bad_inds0, reverse=True):
#     outlier_V0.append(vel0_list.pop(idx))
#     outPh0.append(ph0_list.pop(idx))
#     outTh0.append(th0_list.pop(idx))
# for idx in sorted(bad_inds1, reverse=True):
#     outlier_V1.append(vel1_list.pop(idx))
#     outPh1.append(ph1_list.pop(idx))
#     outTh1.append(th1_list.pop(idx))
# for idx in sorted(bad_inds1b, reverse=True):
#     outlier_V1.append(vel1_list.pop(idx))
#     outPh1.append(ph1_list.pop(idx))
#     outTh1.append(th1_list.pop(idx))

# N_outliers0 = len(outlier_V0)
# N_outliers1 = len(outlier_V1)


# # Convert vel0 and vel1 back to arrays
# vel0 = np.array(vel0_list)
# vel1 = np.array(vel1_list)
# ph0 = np.array(ph0_list)
# th0 = np.array(th0_list)
# ph1 = np.array(ph1_list)
# th1 = np.array(th1_list)
