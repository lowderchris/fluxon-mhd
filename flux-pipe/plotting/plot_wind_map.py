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

from magnetoget import load_magnetogram_params


# create the argument parser
parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
parser.add_argument('--cr', type=int, default=0, help='Carrington Rotation')
parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
parser.add_argument('--show', type=int, default=1)
parser.add_argument('--interp', type=str, default="cubic")
parser.add_argument('--nact', type=int, default=0)
parser.add_argument('--batch', type=str, default="fluxon r2", help='select the batch name')
args = parser.parse_args()

batch = args.batch

(hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)
CR = args.cr or cr or 2183

print(f"plotting CR{CR}...", end="\n" if __name__=="__main__" else "")
filename = f'{args.dat_dir}/{batch}/cr{CR}/wind/radial_wind.dat'

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

# Remove outliers from the dataset
def remove_outliers(data, ph0, th0, threshold=3):
    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    filtered_data, good_points = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data) if mean - threshold * std < x < mean + threshold * std]))
    ph0c, th0c = ph0[good_points], th0[good_points]
    bad_points = [i for i in range(len(data)) if i not in good_points]
    ph0b, th0b = ph0[bad_points], th0[bad_points]
    outlier_data = data[bad_points]
    mean = np.mean(filtered_data)
    std = np.std(filtered_data)
    return np.asarray(filtered_data), mean, std, ph0c, th0c, outlier_data, ph0b, th0b

# Clean the Data
vel0_clean, mean0, std0, ph0_clean, th0_clean, outlier_V0, ph0b, th0b = remove_outliers(vel0, ph0, th0)
vel1_clean, mean1, std1, ph1_clean, th1_clean, outlier_V1, ph1b, th1b = remove_outliers(vel1, ph1, th1)
vel1_clean, mean1, std1, ph1_clean, th1_clean, outlier_V1b, ph1bb, th1bb = remove_outliers(vel1_clean, ph1_clean, th1_clean)

np.asarray(list(outlier_V1).extend(outlier_V1b))
np.asarray(list(ph1b).extend(ph1bb))
np.asarray(list(th1b).extend(th1bb))
# ph1b.extend(ph1bb); 
# th1b.extend(th1bb)


# Plot a histogram of vel0 and vel1 as two subplots
fig, ax = plt.subplots(2, 1, figsize=(8, 8))

ax[0].hist(vel0_clean, bins=10)
ax[0].axvline(mean0, color='k', linestyle='dashed', linewidth=1, label="mean: {:.2f}".format(mean0))
ax[0].axvline(mean0+std0, color='k', linestyle=':', linewidth=1, label="std: {:.2f}".format(std0))
ax[0].axvline(mean0-std0, color='k', linestyle=':', linewidth=1)
ax[0].legend()

ax[1].hist(vel1_clean, bins=10, color='sandybrown')
ax[1].axvline(mean1, color='k', linestyle='dashed', linewidth=1, label="mean: {:.2f}".format(mean1))
ax[1].axvline(mean1+std1, color='k', linestyle=':', linewidth=1, label="std: {:.2f}".format(std1))
ax[1].axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
ax[1].legend()

ax[0].set_title("Near Surface")
ax[1].set_title("Outer Boundary")

ax[0].set_xlabel("Velocity (km/s)")
ax[1].set_xlabel("Velocity (km/s)")

ax[0].set_ylabel("Number of Fluxons")
ax[1].set_ylabel("Number of Fluxons")

fig.suptitle(f"Solar Wind Velocity Distribution: CR{CR}, Open Fluxons: {len(vel1)}", fontsize=16)

plt.tight_layout()
plt.savefig(f"{args.dat_dir}/{batch}/cr{CR}/wind/vel_hist_{len(vel1)}.png")
# plt.show()
plt.close(fig)



# Find the positive max, and min
vel0_max = np.nanmax(vel0_clean) #np.nanmax(np.abs(vel0)) 
vel1_max = np.nanmax(vel1_clean) #np.nanmax(np.abs(vel1)) 
vel0_min = np.nanmin(vel0_clean) #np.nanmin(np.abs(vel0)) 
vel1_min = np.nanmin(vel1_clean) #np.nanmin(np.abs(vel1)) 


scale = 15**2
power = 1

v0 = scale * ((np.abs(vel0_clean) - vel0_min) / (vel0_max-vel0_min))**power
v1 = scale * ((np.abs(vel1_clean) - vel1_min) / (vel1_max-vel1_min))**power

outlier_V0_scaled = scale * ((np.abs(outlier_V0) - vel0_min) / (vel0_max-vel0_min))**power
outlier_V1_scaled = scale * ((np.abs(outlier_V1) - vel1_min) / (vel1_max-vel1_min))**power


which_interp = ['nearest', 'cubic', 'linear'] if args.interp else ['linear']
for interp in which_interp:
    # Plot the Data
    fig, (ax0, ax1, ax2) = plt.subplots(3, sharex = True, sharey= True )

    magnet = load_fits_magnetogram()

    # find the max and min of the magnetogram plot for use in setting the colormap, 
    sigma = 3
    mmean = np.nanmean(magnet)
    msig = np.nanstd(magnet)
    mvmin = mmean - sigma*msig
    mvmax = mmean + sigma*msig


    # Plot the magnetogram
    magnetic_vminmax = [-500, 500]
    magimg = ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), 
                        aspect='auto', vmin=mvmin, vmax=mvmax)
    # Scatter the fluxon data
    ph000 = (ph0_clean+np.pi)%(2*np.pi)
    ph111 = (ph1_clean+np.pi)%(2*np.pi)
    ph000b = (ph0b+np.pi)%(2*np.pi)
    ph111b = (ph1b+np.pi)%(2*np.pi)

    sc00 = ax0.scatter(ph000, th0_clean, c= vel0_clean, s=v0, alpha=0.5, label='V(1.0R$_\odot$)'  , cmap="winter",  )
    sc01 = ax0.scatter(ph111, th1_clean, c= vel1_clean, s=v1, alpha=0.5, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
    # sc01 = ax1.scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
    # sc10 = ax1.scatter(ph0, th0, c= fr0 , s=f0, alpha=0.75, label='Fr(1.0R$_\odot$)' , cmap="winter", )
    # sc11 = ax1.scatter(ph1, th1, c= fr1 , s=f1, alpha=0.75, label='Fr(21.5R$_\odot$)', cmap="autumn", marker='s')

    #Scatter the Outliers
    ax0.scatter(ph000b, th0b, c= outlier_V0_scaled, alpha=0.6, label='V(1.0R$_\odot$)', cmap="winter", marker='X',  s=50, edgecolors='k')
    ax0.scatter(ph111b, th1b, c= outlier_V1_scaled, alpha=0.6, label='V(1.5R$_\odot$)', cmap="autumn", marker='X',  s=50, edgecolors='k')
    ax1.scatter(ph000b, th0b, c= outlier_V0_scaled, alpha=0.95, label='V(1.0R$_\odot$)', cmap="winter", marker='X', s=50, edgecolors='k')
    ax2.scatter(ph111b, th1b, c= outlier_V1_scaled, alpha=0.95, label='V(1.5R$_\odot$)', cmap="autumn", marker='X', s=50, edgecolors='k')



    from scipy.interpolate import griddata

    ## Plot the interpolated data
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
    sc00 = ax1.scatter(ph0, th0, c=vel0, s=v0, alpha=1, label='V(1.0R$_\odot$)', cmap="winter", edgecolors='k'             )
    sc01 = ax2.scatter(ph1, th1, c=vel1, s=v1, alpha=1, label='V(21.5R$_\odot$)', cmap="autumn", marker='s', edgecolors='k')

    # Create a contour plot
    contour0 = ax1.contourf(grid_x, grid_y, grid_z0, zorder=0, alpha=1, cmap="winter")
    contour1 = ax2.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap="autumn")



    # Add a colorbar
    cbar0 = fig.colorbar(contour0, ax=ax1)
    cbar1 = fig.colorbar(contour1, ax=ax2)
    cbar0.set_label('Interpolated Wind [km/s]\nV( 1.0R$_\odot$)')
    cbar1.set_label('Interpolated Wind [km/s]\nV(21.5R$_\odot$)')

    # ax0.set_


    cbar01 = fig.colorbar(magimg, ax=ax0, extend='both')
    cbar01.set_label("Mag Field Anomaly [Gauss]")

    # cbar01 = fig.colorbar(sc01, ax=ax0)
    # cbar00 = fig.colorbar(sc00, ax=ax0)
    # cbar11 = fig.colorbar(sc11, ax=ax1)
    # cbar10 = fig.colorbar(sc10, ax=ax1)

    # cbar00.set_label("V(1.0R$_\odot$)  [km/s]")
    # cbar01.set_label("V(21.5R$_\odot$) [km/s]")
    # cbar10.set_label("Fr(1.0R$_\odot$)  ")
    # cbar11.set_label("Fr(21.5R$_\odot$) ")


    fig.suptitle(F"Fluxon Wind Strength for CR {CR}: N={len(v1)}")
    # ax1.set_title(F"Fluxon Expansion for CR {args.cr}")
    for ax in [ax0, ax1, ax2]:
        ax.set_ylabel('Sine latitude')
        ax.set_ylim((-1.5,1.1))
        ax.axhline(-1, c='lightgrey', zorder=-10)
        ax.axhline( 1, c='lightgrey', zorder=-10)
        ax.axvline(0, c='lightgrey', zorder=-10)
        ax.axvline(2*np.pi, c='lightgrey', zorder=-10)

    ax2.set_xlabel('Longitude (Radians)')
    # ax0.legend(loc="lower right")
    ax0.set_ylim((-1, 1))
    ax1.set_ylim((-1, 1))
    ax2.set_ylim((-1, 1))

    fig.set_size_inches((8,8))
    plt.tight_layout()
    pngname = filename.replace(".dat", f"_map_{len(vel1)}_{interp}.png")
    plt.savefig(pngname, dpi=200)

    import os.path as path
    import os

    file = f"cr{CR}_radial_wind_map_ou{len(vel1)}_f{args.nact}.png"
    directory = f"{args.dat_dir}/{batch}/imgs/wind_maps"
    if not path.exists(directory):
        os.makedirs(directory)
        
    save_path = path.join(directory, file)
    # print("Saving fluxon image to ", save_path)
    print(save_path)
    plt.savefig(save_path, dpi=150)

    # if args.show:
    #     plt.show()
    plt.close(fig)


    def plot_three(grid_x, grid_y, grid_z, save_path):
        fig, axarr = plt.subplots(3, sharex="all", sharey="all", figsize=(10,10))
        # the imshow extent will let it plot on the same space as the other two plots, using grid_x and grid_y
        extent = (np.min(grid_x), np.max(grid_x), np.min(grid_y), np.max(grid_y))
        axarr[0].set_title("imshow")
        axarr[0].imshow(grid_z, cmap ="autumn", extent=extent, interpolation='None', origin='lower')
        axarr[1].set_title("contourf")
        axarr[1].contourf(grid_x, grid_y, grid_z, zorder=0, alpha=1, cmap="autumn")
        axarr[2].set_title("hexbin")

        nx = 12
        sy, sx = grid_x.shape
        ratio = float(sy)/float(sx)
        gridsize = (nx, int(np.round(nx*ratio)))

        axarr[2].hexbin(grid_x.flatten(), grid_y.flatten(), C=grid_z.flatten(), gridsize=gridsize, cmap='autumn')
        
        for ax in axarr:
            ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig)
        # plt.show()

    pngname = filename.replace(".dat", f"_map_{len(vel1)}_{interp}_options.png")

    # plot_three(grid_x, grid_y, grid_z1, pngname)

# ax0.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*vel0/vel0.max())**3, alpha=0.75, label='V(1.0R$_\odot$)')
# ax0.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*vel1/vel1.max())**3, alpha=0.75, label='V(21.5R$_\odot$)')
# ax1.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*fr0/fr0.max())**3, alpha=0.75, label='Fr(1.0R$_\odot$)', marker='s')
# ax1.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*fr1/fr1.max())**3, alpha=0.75, label='Fr(21.5R$_\odot$)', marker='s')

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
