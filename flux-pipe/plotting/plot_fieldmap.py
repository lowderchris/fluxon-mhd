

# datdir = "/Users/cgilbert/vscode/fluxon-data/"



# read_fr.py
# Script to parse fluxon output file of the field expansion factor


# file output format: fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1
# fluxon id, beginning coords, end coords, beginning mag, end mag, area maybe?




import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
import os.path as path
import os
import sys
import plot_helper

from pfss_funcs import pixel_to_latlon
from magnetoget import load_fits_magnetogram, load_magnetogram_params
from magnetoget import shorten_path

# print("\n \n \n Plotting FIELD MAP...", end="")
# create the argument parser
def get_ax(ax=None):
    if ax is not None:
        fig, ax0 = ax.get_figure(), ax
    else:
        fig, ax0 = plt.subplots(1)
    return fig, ax0

def magnet_plot(get_cr, datdir, batch, open_f=None, closed_f=None, force=False, reduce=0, nact=0, nwant=None, do_print_top=False,
                ax=None, verb=True, ext = "pdf", plot_all=True, plot_open=True, do_print=False, vmin=-500, vmax=500):
    
    fig, ax0 = get_ax(ax)
    # if do_print:
    #     print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #     print("(py) Plotting Magnetic Field Maps")
    #     print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


    if do_print_top: print("\n\t\tMaking Magnetogram with Footpoints...")
    # Define the directory paths for the files
    floc_path = f"{datdir}/batches/{batch}/cr{get_cr}/floc/"
    top_dir   = f"{datdir}/batches/{batch}/imgs/footpoints/"
    if not path.exists(top_dir):
        os.makedirs(top_dir)

    # Define the file names with their complete paths
    open_file   = open_f     or   f"{floc_path}floc_open_cr{get_cr}_r{reduce}_f{nwant}.dat"
    closed_file = closed_f or     f"{floc_path}floc_closed_cr{get_cr}_r{reduce}_f{nwant}.dat"
    all_file    = closed_file.replace("closed_", "")
    # print(open_f)
    # print(closed_f)
    # Load the data
    if do_print_top: print(f"\t\t\tOpening {shorten_path(all_file)}...")
    fluxon_location = np.genfromtxt(all_file)
    
    if do_print_top: print(f"\t\t\tOpening {shorten_path(open_file)}...")
    oflnum, oflx, olat, olon, orad = np.loadtxt(open_file, unpack=True)

    if do_print_top: print(f"\t\t\tOpening {shorten_path(closed_file)}...")
    cflnum, cflx, clat, clon, crad = np.loadtxt(closed_file, unpack=True)


    from magnetoget import load_fits_magnetogram, load_magnetogram_params

    magnet, header = load_fits_magnetogram(batch=batch, ret_all=True)
    f_lat, f_lon, f_sgn, fnum = pixel_to_latlon(magnet, header, fluxon_location)
        # # ph0, th0 = phi0+np.pi, -(theta0-(np.pi/2))
        # # ph1, th1 = phi1+np.pi, -(theta1-(np.pi/2))

    ## Keep only the values where the radius is 1.0
    rtol = 0.01
    get_r = 1.0
    #Open fields
    oflnum_low = oflnum[np.isclose(orad, get_r, rtol)]
    oflx_low =     oflx[np.isclose(orad, get_r, rtol)]
    olat_low =     olat[np.isclose(orad, get_r, rtol)]
    olon_low =     olon[np.isclose(orad, get_r, rtol)]
    orad_low =     orad[np.isclose(orad, get_r, rtol)]

    # Closed fields
    cflnum_low = cflnum[np.isclose(crad, get_r, rtol)]
    cflx_low =     cflx[np.isclose(crad, get_r, rtol)]
    clat_low =     clat[np.isclose(crad, get_r, rtol)]
    clon_low =     clon[np.isclose(crad, get_r, rtol)]
    crad_low =     crad[np.isclose(crad, get_r, rtol)]

    # Convert to radians
    ph_olow, th_olow = np.sin(np.deg2rad(olat_low)), np.deg2rad(olon_low)
    ph_clow, th_clow = np.sin(np.deg2rad(clat_low)), np.deg2rad(clon_low)

    # Report the number of open and closed fluxons
    n_open = int(np.max(oflnum_low))
    n_closed = int(np.max(cflnum_low))
    n_flux = n_open + n_closed
    n_outliers = np.abs(fnum-n_flux)

    
    # Define the file name for the plot
    pic_name = f'cr{get_cr}_f{fnum}_ou{n_open}_footpoints_topology.{ext}'
    fluxon_map_output_path =   path.join(floc_path, pic_name)
    fluxon_map_output_path_top = path.join(top_dir, pic_name)

    # Check if the plot already exists
    do_plot = False
    pic_paths = [fluxon_map_output_path, fluxon_map_output_path_top]
    for testpath in pic_paths:
        if not path.exists(testpath):
            do_plot = True
            break

    if do_plot or force or (ax is not None):
        if do_print: print(f"\tPlotting...", end="")
        # Plot ###################
        # Plot the magnetogram
        sigma = 2
        mmean = np.nanmean(magnet)
        msig = np.nanstd(magnet)
        mvmin = mmean - sigma*msig
        mvmax = mmean + sigma*msig
        ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto', 
                #    vmin=mvmin, vmax=mvmax)
                   vmin=vmin, vmax=vmax)

        # Plot all the fluxons
        if plot_all:
            colors_all = ['orange' if s > 0 else 'teal' for s in f_sgn]
            ax0.scatter(f_lon, f_lat, s=3**2, c=colors_all, alpha=0.6)

        # Plot the open fluxons
        if plot_open:
            colors_open = ['red' if s > 0 else 'blue' for s in oflx_low]
            # sc00 = ax0.scatter(np.deg2rad(olon_low), np.sin(olat_low), c=oflx_low, cmap="PiYG", alpha=0.75, label='B(1.0R$_\odot$)')
            # sc00 = ax0.scatter(th_clow, ph_clow, c=cflx_low, cmap="PuOr_r", alpha=0.55, label='B(1.0R$_\odot$)')
            sc00 = ax0.scatter(th_olow, ph_olow, s=5**2, c=colors_open, alpha=1.0, label='Open Fields', edgecolors='k')

        if ax is None:
            shp = magnet.shape #pixels
            plt.axis('off')
            sz0=6 #inches
            ratio = shp[1]/shp[0]
            sz1=sz0*ratio #inches
            DPI = shp[1] / sz1 #pixels/inch
            fig.set_size_inches((sz1, sz0))
            plt.tight_layout()
            print(f"\tSaving {fluxon_map_output_path}...")
            plt.savefig(fluxon_map_output_path, bbox_inches='tight', dpi=4*DPI)
            print(f"\tSaving {fluxon_map_output_path_top}...")
            plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)
            # plt.show()
            plt.close(fig)
    else:
        if do_print: 
            print(f"Skipped! Files already exist:")
            print(f"\t    {shorten_path(fluxon_map_output_path)}")
            print(f"\t    {shorten_path(fluxon_map_output_path_top)}")
    if do_print: print(f"\n\tn_open: {n_open}, n_closed: {n_closed}, n_total: {n_flux}, n_all: {fnum}, n_outliers: {n_outliers}\n")
    
    if do_print_top: print("\t\t    Success!")

    return n_open, n_closed, n_flux, fnum, n_outliers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None)
    parser.add_argument('--nwant', type=int, default=None)
    parser.add_argument('--open', type=str, default=None)
    parser.add_argument('--closed', type=str, default=None)
    parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
    parser.add_argument('--batch', type=str, default="default_batch", help='select the batch name')
    args = parser.parse_args()
    batch = args.batch

    (hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)
    CR = args.cr or cr
    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(CR, args.dat_dir, batch, args.open, args.closed, do_print=True, reduce=reduce, nwant=args.nwant)
    
    # exit()

# print("DONE!")
# # sc00 = ax0.scatter(np.deg2rad(clon_low), np.sin(clat_low), c=cflx_low, cmap="PuOr_r", alpha=0.75, label='B(1.0R$_\odot$)')

# # sc01 = ax0.scatter(clat, clon, c=br1, s = b1, cmap="autumn", alpha=0.75, label='B(21.5R$_\odot$)', marker='s')



# filename = f'{args.dat_dir}/{batch}/cr{CR}/wind/radial_bmag.dat'



# # Load the dat file
# arr = np.loadtxt(filename).T
# fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1 = arr
# nfluxon = arr.shape[1]


# # Convert coords to correct coords
# ph0, th0 = phi0+np.pi, -(theta0-(np.pi/2))
# ph1, th1 = phi1+np.pi, -(theta1-(np.pi/2))


# # Do some data manipulation
# br0_max = np.nanmax(br0) or 0.25
# br1_max = np.nanmax(br1) or 0.25
# ar0_max = np.nanmax(ar0) or 0.25
# ar1_max = np.nanmax(ar1) or 0.25

# skew = 5**2
# power = 1
# b0 = skew*(6*br0/br0_max)**power
# b1 = skew*(4*br1/br1_max)**power
# a0 = skew*(6*ar0/ar0_max)**power
# a1 = skew*(4*ar1/ar1_max)**power


# # Plot the Data
# fig, (ax0, ax1) = plt.subplots(2)

# magnet = load_fits_magnetogram()
# ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')
# ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')

# sc00 = ax0.scatter(ph0, th0, c=br0, s = b0, cmap="winter", alpha=0.75, label='B(1.0R$_\odot$)')
# sc01 = ax0.scatter(ph1, th1, c=br1, s = b1, cmap="autumn", alpha=0.75, label='B(21.5R$_\odot$)', marker='s')
# sc10 = ax1.scatter(ph0, th0, c=ar0, s = a0, cmap="winter", alpha=0.75, label='A(1.0R$_\odot$)')
# sc11 = ax1.scatter(ph1, th1, c=ar1, s = a1, cmap="autumn", alpha=0.75, label='A(21.5R$_\odot$)', marker='s')

# cbar01 = fig.colorbar(sc01, ax=ax0)
# cbar00 = fig.colorbar(sc00, ax=ax0)
# cbar11 = fig.colorbar(sc11, ax=ax1)
# cbar10 = fig.colorbar(sc10, ax=ax1)

# cbar00.set_label("B(1.0R$_\odot$)  [Gauss] ")
# cbar01.set_label("B(21.5R$_\odot$) [Gauss]")
# cbar10.set_label("A(1.0R$_\odot$)  [m$^2$]")
# cbar11.set_label("A(21.5R$_\odot$) [m$^2$]")

# ax0.set_title(F"Fluxon Magnetic Field Strength for CR {get_cr}")
# ax1.set_title(F"Fluxon Area for CR {get_cr}")

# for ax in (ax0, ax1):
#     ax.set_xlabel('Longitude (Radians)')
#     ax.set_ylabel('Sine latitude')
#     ax.set_ylim((-1.5,1.1))
#     ax.axhline(-1, c='lightgrey', zorder=-10)
#     ax.axhline( 1, c='lightgrey', zorder=-10)
#     ax.axvline(0, c='lightgrey', zorder=-10)
#     ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
#     ax.legend(loc="lower right")

# fig.set_size_inches((8,8))
# plt.tight_layout()
# pngname = filename.replace(".dat", f"_{len(ph0)}_{len(ph1)}.png")
# # pngname = filename.replace(".dat", ".png")
# plt.savefig(pngname)
# if args.show:
#     plt.show()
# plt.close(fig)
# print("Done!")





# ax0.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*br0/br0_max)**3, alpha=0.75, label='B(1.0R$_\odot$)')
# ax0.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*br1/br1_max)**3, alpha=0.75, label='B(21.5R$_\odot$)')
# ax1.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*ar0/ar0_max)**3, alpha=0.75, label='A(1.0R$_\odot$)', marker='s')
# ax1.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*ar1/ar1_max)**3, alpha=0.75, label='A(21.5R$_\odot$)', marker='s')