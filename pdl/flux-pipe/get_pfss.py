"""
Fluxon PFSS mapping
=======================
Generates a fluxon mapping from input GONG-sourced pfss coronal field solution
"""


datdir = "/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/"


###############################################################################
# First, import required modules
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import os
import sys
import astropy.constants as const
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sunpy.map
from scipy import ndimage, misc
import scipy.interpolate
import pfsspy
from pfsspy import coords
from pfsspy import tracing
import pickle
from magnetoget import load_magnetogram_params, read_fits_data

# Helper Function to Read the Fits Files

if len(sys.argv)>1: 
    # print(int(sys.argv[1])>0)
    do_pfss = force_plot = force_trace = True if int(sys.argv[1])>0 else False
else:
    do_pfss = force_plot = force_trace = False

# print(do_pfss, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

# Load the parameters, basically just which magnetogram is used
print("\n -->Loading Parameters...")
(hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(datdir)


###############################################################################
# We can now use SunPy to load the HMI fits file, and extract the magnetic
# field data. Interpolate the data down to a more reasonable resolution.
#
# The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is
# not obvious this is the correct way to do this, so use the following lines
# at your own risk!

# brmap = (sunpy.io.fits.read(datdir + fname)) # b radial map
# brdat = brmap[0].data[2,:,:] # b data itself

# Load the fits file and format the data and header
print("\n -->Loading and conditioning fits file...")
fits_path = datdir + fname
hdulist = read_fits_data(fits_path)
brdat = hdulist[0].data
header= hdulist[0].header

brdat = brdat - np.mean(brdat)
br = sunpy.map.Map(brdat, header)

br.meta['ctype1'] = 'CRLN-CEA'
br.meta['ctype2'] = 'CRLT-CEA'
br.meta['naxis'] = 2

br.meta['crlt_obs'] = 0.0
br.meta['crln_obs'] = 0.0
br.meta['dsun_obs'] = 1.0
# br.meta['bunit'] = "Gauss"
# if "rsun_obs" not in br.meta:
#     from astropy import units as u
#     import sunpy.map
#     from sunpy.data.sample import AIA_171_IMAGE
#     aiamap = sunpy.map.Map(AIA_171_IMAGE)
#     # br.meta['rsun_obs'] = (696340 * u.km / aiamap.dsun).to(u.arcsec).value
#     # Set the solar radius metadata (using a value in arcseconds)
#     br.meta["rsun_obs"] = ((696340 * u.km))

# br.meta['cdelt1'] = 360 / brdat.shape[1]
# br.meta['cdelt2'] *= int(reduce)
# br.meta['cdelt2'] = 360 / (brdat.shape[0] * np.pi)

import copy
br_safe = copy.copy(br)
br_safe.data[np.isnan(br_safe.data)]=0 

###############################################################################
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.
nrho = 60
rss = 2.5

###############################################################################
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an Input object that stores this information
print("\n -->Getting Pfss...", end="")
import os.path as path
pickle_dir = path.dirname(path.dirname(fits_path))
pickle_path = path.join(pickle_dir, "pfss_output.pkl")

load_failure = False
if not do_pfss:
    print("Loading from Disk...", end="")
    # Load the results
    try:
        with open(pickle_path, 'rb') as inp:
            output = pickle.load(inp)
        print("Success!")
    except FileNotFoundError as e:
        print("File not found!")
        # print(e, file=sys.stderr,)
        load_failure = True
        print("\n -->Computing PFSS...", end="")
        
###############################################################################
# Now calculate the PFSS solution
if do_pfss or load_failure:
    input = pfsspy.Input(br_safe, nrho, rss)

    # An attempt to use the performance cores
    # from os import setpriority, getpid, getpriority
    # PRIO_DARWIN_THREAD  = 0b0011
    # PRIO_DARWIN_PROCESS = 0b0100
    # PRIO_DARWIN_BG      = 0x1000
    # setpriority(PRIO_DARWIN_PROCESS, 0, 0)

    # Actually do the math
    output = pfsspy.pfss(input)
    print("Done! Saving Pickled Result.")

    # Save the results
    with open(pickle_path, 'wb') as outp:
        pickle.dump(output, outp, pickle.HIGHEST_PROTOCOL)

    
# output.bg
# output.bc
# output.bunit

#bc: B on the centres of the cell faces.
#bg: B as a (weighted) averaged on grid points.
#bunit: Unit of the input map data.
###############################################################################


# CL - DOUBLE CHECK THE LATITUTE DIRECTION HERE, and plot the polarity inversion line.
fluxon_location = np.genfromtxt(datdir + 'fluxon/cr' + cr + '/floc/floc_cr'+cr+'.dat')
shp = br_safe.data.shape
n_lat = shp[0]
lat_center = n_lat // 2
f_lon = np.deg2rad(fluxon_location[:,0] * br_safe.meta['cdelt1'])
f_lat = (fluxon_location[:,1]-lat_center) * br_safe.meta['cdelt2']
f_sgn = fluxon_location[:,2]
n_flux = len(f_sgn)

fluxon_map_output_path = path.join(path.dirname(path.dirname(fits_path)), f'{n_flux}_footprint.png')

print("\n -->Plotting Fluxon Locs...", end="")
if not path.exists(fluxon_map_output_path) or force_plot:
    # br = sunpy.map.Map(datdir + 'hmi.Synoptic_Mr.polfil/hmi.synoptic_mr_polfil_720s.' + cr + '.Mr_polfil.fits')
    # fluxon_map_output_path = fits_path.replace('.fits', '_fluxons.png')
    # fluxon_map_output_path = path.join(path.dirname(path.dirname(fits_path)), path.basename(fits_path).replace('.fits', f'_{n_flux}_fluxons.png'))

    ## Print the Fluxon Map
    fig, ax = plt.subplots()
    ax.imshow(br_safe.data, cmap='gray', interpolation=None, origin="lower", zorder=-10)
    ratio = shp[1]/shp[0]
    for (x, y, sig) in fluxon_location:
        color = 'red' if sig > 0 else "teal"
        ax.scatter(x, y, c=color, alpha=0.4)

    # import pdb
    # pdb.set_trace()
    # ax.imshow(output.bc[0], cmap='RdBu', interpolation=None, origin="lower", zorder=-5, alpha= 0.33)
    # pfss_out = output
    # ss_br = pfss_out.source_surface_br
    # # Create the figure and axes
    # # fig2 = plt.figure()
    # ax = plt.subplot(projection=ss_br)
    # fig, ax = plt.subplots()
    # # Plot the source surface map
    # ss_br.plot()
    # # Plot the polarity inversion line
    # ax.plot_coord(pfss_out.source_surface_pils[0])
    # # Plot formatting
    # plt.colorbar()
    # ax.set_title('Source surface magnetic field')

    # plt.show()
    # output.bc
    # output.bg
    # output.bunit

    plt.axis('off')
    sz0=6 #inches
    sz1=sz0*ratio #inches
    shp #pixels
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.tight_layout()
    plt.savefig(fluxon_map_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()

    ## Plot the fluxon map in lat/lon
 
    fig, ax = plt.subplots()
    ax.imshow(br_safe.data, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')

    for (long, lat, sig) in zip(f_lon, f_lat, f_sgn):
        color = 'red' if sig > 0 else "teal"

        ax.scatter(long, lat, color=color)
    # fluxon_map_output_path2 = path.join(path.dirname(path.dirname(fits_path)), path.basename(fits_path).replace('.fits', f'_latlon_{n_flux}_fluxons.png'))
    fluxon_map_output_path2 = path.join(path.dirname(path.dirname(fits_path)), f'{n_flux}_footprints_latlon.png')
    plt.axis('off')
    fig.set_size_inches((sz1, sz0))
    plt.tight_layout()
    plt.savefig(fluxon_map_output_path2, bbox_inches='tight', dpi=4*DPI)
    plt.close()
    # plt.show()
    print("Success!")
else:
    print("Skipped! File Already Exists")

print("\n -->Tracing Open and Closed Fluxons...")

# Note that this code was originally developed using an older version of the pfsspy code - improvements look to be able to handle bundles of fieldlines, which would greatly simplify some of the code. Do this.

open_path = datdir + 'fluxon/cr' + cr + '/floc/' + 'floc_open_cr'+cr+'.dat'
closed_path = datdir + 'fluxon/cr' + cr + '/floc/' + 'floc_closed_cr'+cr+'.dat'


from tqdm import tqdm
import timeout_decorator


# Capture an array of starting and endpoints for open fieldlines 2x(r, lat, lon) + sgn of starting point
fl_open = np.zeros([1,5])
flnum_open = 0

# Capture an array of starting and endpoints for closed fieldlines 2x(r, lat, lon) + sgn of starting point
fl_closed = np.zeros([1,5])
flnum_closed = 0

r0 = 1.01 * const.R_sun
coord_frame = output.coordinate_frame

@timeout_decorator.timeout(1)
def trace_each(coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer):
    #x0 = np.array(coords.sph2cart(r0, f_lat[i], flon[i]))
    (this_flon, this_flat) = coords
    x0 = SkyCoord(this_flon * u.rad, this_flat * u.rad, r0, frame=coord_frame)
    fl = output.trace(tracer, x0)
    fl = fl.field_lines[0]
    if fl.is_open:
        fl_pol = fl.polarity
        fl = fl.coords.spherical
        fl_lats = fl.lat.value
        fl_lons = fl.lon.value
        fl_rads = (fl.distance/const.R_sun).value
        if (fl_pol == -1): (fl_lats, fl_lons, fl_rads) = (fl_lats[::-1], fl_lons[::-1], fl_rads[::-1])
        fl_lats = np.concatenate((fl_lats, fl_lats[-1]*np.ones(10)),axis=0)
        fl_lons = np.concatenate((fl_lons, fl_lons[-1]*np.ones(10)),axis=0)
        fl_rads = np.concatenate((fl_rads, np.linspace(2.5,22,10)),axis=0) 
        prev_rad = 0.
        for j in np.arange(0, len(fl_lats)):
            # Output flnum, polarity, latitude, longitude, radius
            if np.abs(fl_rads[j] - prev_rad) > 0.1:
                fl_open = np.append(fl_open, [[flnum_open, f_sgn[i], fl_lats[j], fl_lons[j], fl_rads[j]]], axis=0)
                prev_rad = fl_rads[j]
        flnum_open += 1
    else:
        fl = fl.coords.spherical
        fl_rads = (fl.distance/const.R_sun).value
        max_rad = fl_rads.max()
        rad_thresh = np.max([((max_rad-1.) / 5.), 0.01])
        prev_rad = 0.
        if np.mod(flnum_closed, 2) == 0:
            for j in np.arange(0, len(fl.lat)):
                # Output flnum, polarity, latitude, longitude, radius
                if (np.abs(fl_rads[j] - prev_rad) > rad_thresh) or (fl_rads[j] == max_rad):
                    fl_closed = np.append(fl_closed, [[flnum_closed, f_sgn[i], fl.lat[j].value, fl.lon[j].value, (fl.distance[j]/const.R_sun).value]], axis=0)
                    prev_rad = fl_rads[j]
        flnum_closed += 1
    return output, fl_open, fl_closed, flnum_open, flnum_closed, tracer

def trace_lines(output, f_lon, f_lat, fl_open, fl_closed, flnum_open, flnum_closed):
    skip_num = 0
    timeout_num = 0
    tracer = tracing.PythonTracer()
    for i, coords in enumerate(tqdm(zip(f_lon, f_lat), desc="Tracing Field Lines", total=len(f_lat))):
        try:
            output, fl_open, fl_closed, flnum_open, flnum_closed, tracer = trace_each(coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer)
        except timeout_decorator.TimeoutError as e:
            # print(f"Timeout Occured in Iteration {i}")
            timeout_num += 1
        except ValueError as e:
            # print(f"Skipping {i}: ",e)
            # try:
            #     (lon, lat) = coords
            #     new_lat = lat if lat>0 else lat
            #     new_coords = (lon, new_lat)
            #     fl_open, fl_closed = trace_each(new_coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed)
            # except ValueError as e:
            skip_num += 1

    if skip_num > 0 or timeout_num > 0:
        t_perc = 100*timeout_num/len(f_lon)
        s_perc = 100*skip_num/len(f_lon)
        print(f"\n\nSome iterations failed. Timed-out: {timeout_num} ({t_perc:0.2f}%), ValueError: {skip_num} ({s_perc:0.2f}%)\n\n")

    print(f"\n\tOpen Lines: {len(fl_open)}, Closed Lines: {len(fl_closed)}, Failures: {skip_num+timeout_num}")
        
    fl_open = fl_open[1:]
    fl_closed = fl_closed[1:]
    return fl_open, fl_closed
    
if not os.path.exists(open_path) or force_trace:
    fl_open, fl_closed = trace_lines(output, f_lon, f_lat, fl_open, fl_closed, flnum_open, flnum_closed)

    # Output is flnum, polarity, latitude, longitude, radius

    # Save these coordinates to file
    print("\n -->Saving Fluxons...", end="")
    np.savetxt(open_path, fl_open)
    np.savetxt(closed_path, fl_closed)
    print("Success!")
else:
    print("\tSkipped! floc dat files already exist.")





# # print("\n-->>>>>>>>>>>>>>\n-->>>>Main Program Complete!<<<<\n<<<<<<<<<<<<<<\n")
# print(" -->Plotting: ", bool(doplot), end="\n\n")
# def set_axes_lims(ax):
#     ax.set_xlim(0, 360)
#     ax.set_ylim(0, 180)

# # If you want to plot...
# doplot = False
# if doplot and False:
#     print("\n -->Plotting...")

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     r = 1.01
#     for i in np.arange(0,len(f_lon)):
#         x0 = SkyCoord(f_lon[i] * u.rad, f_lat[i] * u.rad, r)

#         # x0 = coords.sph2cart(r, f_lat[i], f_lon[i])
#         fl = output.trace(tracer, x0)
#         fl = fl.field_lines[0]
#         color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(fl.polarity)
#         ax.plot(fl.x / const.R_sun,
#                 fl.y / const.R_sun,
#                 fl.z / const.R_sun,
#                 color=color, linewidth=1)
#         if fl.is_open:
#             fl = fl.spherical
#             x0 = np.array(coords.sph2cart(2.45, -1*(-1*fl.lat[-1].value-90)*np.pi/180, fl.lon[-1].value*np.pi/180))
#             fl = output.trace(x0)
#             color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(fl.polarity)
#             ax.plot(fl.x / const.R_sun,
#                     fl.y / const.R_sun,
#                     fl.z / const.R_sun,
#                     color=color, linewidth=1)
#     #ax.set_title('CR'+cr+' PFSS field')
#     plt.show()

# from matplotlib.pyplot import scatter, plot, imshow, xlim, ylim, xlabel, ylabel
# from numpy import where

# if doplot:
#     scatter(fl_open[:,2], fl_open[:,1])
#     scatter(fl_open[:,5], fl_open[:,4])
#     for i in np.arange(0, len(fl_open[:,0])):
#        plot([fl_open[i,2], fl_open[i,5]],[fl_open[i,1], fl_open[i,4]],'k')

#     scatter(fl_closed[:,2], fl_closed[:,1])
#     scatter(fl_closed[:,5], fl_closed[:,4])
#     for i in np.arange(0, len(fl_closed[:,0])):
#        plot([fl_closed[i,2], fl_closed[i,5]],[fl_closed[i,1], fl_closed[i,4]],'k')

#     wp = where(f_sgn>0)[0]
#     wn = where(f_sgn<0)[0]
#     scatter(f_lon[wn], f_lat[wn])
#     scatter(f_lon[wp], f_lat[wp])

#     wp = where(fl_open[:,6]>0)[0]
#     wn = where(fl_open[:,6]<0)[0]
#     scatter(fl_open[wn,2], fl_open[wn,1])
#     scatter(fl_open[wp,2], fl_open[wp,1])

#     imshow(br, vmin=-25, vmax=25, extent=[0,360,-90,90], aspect='auto')
#     scatter(f_lon[where(f_sgn>0)]*180/np.pi, -1*f_lat[where(f_sgn>0)]*180/np.pi+90, c='white')
#     scatter(f_lon[where(f_sgn<0)]*180/np.pi, -1*f_lat[where(f_sgn<0)]*180/np.pi+90, c='black')
#     xlim([0,360])
#     ylim([-90,90])
#     xlabel('Carrington longitude')
#     ylabel('Latitude')

#     cmap = plt.get_cmap('tab10')
#     imshow(br, vmin=-20, vmax=20, extent=[0,360,-1,1], aspect='auto')
#     scatter(fl_closed[:,2], np.sin(fl_closed[:,1]*np.pi/180), s=5.0, color=cmap(0))
#     scatter(fl_closed[:,5], np.sin(fl_closed[:,4]*np.pi/180), s=5.0, color=cmap(1))
#     scatter(fl_open[:,2], np.sin(fl_open[:,1]*np.pi/180), s=5.0, color=cmap(3))
#     xlim([0,360])
#     ylim([-1,1])
#     xlabel('Carrington longitude (CR2193)')
#     ylabel('Sine latitude')

#     plt.show()


# if __name__ == "__main__":
#     sys.exec("python3 mag_runner.py")

