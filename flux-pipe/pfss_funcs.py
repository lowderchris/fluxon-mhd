


from magnetoget import read_fits_data
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import pfsspy
import pickle
from time import time
from os import path
import os
import numpy as np
import sunpy.map
import copy
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import sunpy.map
from magnetoget import load_fits_magnetogram, load_magnetogram_params, shorten_path




def load_and_condition_fits_file(fname, datdir):
    print("\n\tLoading and conditioning fits file...")
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

    try:
        fits_path = fname
        # print(fname)
        hdulist = read_fits_data(fits_path)
    except FileNotFoundError as e:
        fits_path = os.path.join(datdir, fname)
        hdulist = read_fits_data(fits_path)
    print("\t\t", shorten_path(fits_path, 2))


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

    br_safe = copy.copy(br)
    br_safe.data[np.isnan(br_safe.data)]=0

    return br_safe, fits_path


def fits_path_to_pickle_path(fits_path, reduce):
    the_dir = path.dirname(path.dirname(fits_path))
    pickle_path = path.join(the_dir, f"pfss_output_{reduce}.pkl")
    return pickle_path


def load_pfss(pickle_path):
    print("\n\tGetting Pfss...", end="")

    try:
        with open(pickle_path, 'rb') as inp:
            output = pickle.load(inp)
            print(f"Success! Loaded: \n\t\t {shorten_path(pickle_path, 2)}")
        return output
    except FileNotFoundError as e:
        print("File not found.")
        return None


def compute_pfss(br_safe, pickle_path, nrho=60, rss=2.5):
    ###############################################################################
    # The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
    # rho = ln(r), and r is the standard spherical radial coordinate. We need to
    # define the number of rho grid points, and the source surface radius.
    elapsed = 0

    shp = br_safe.data.shape
    print(f"\t\tComputing PFSS on {shp} magnetogram...", end="", flush=True)
    before = time()
    
    ###############################################################################
    # From the boundary condition, number of radial grid points, and source
    # surface, we now construct an Input object that stores this information

    input = pfsspy.Input(br_safe, nrho, rss)

    ###############################################################################
    # Now calculate the PFSS solution
    output = pfsspy.pfss(input)
    elapsed = time() - before
    print("Done! Took {:.2f} seconds.".format(elapsed))

    print("\t\tSaving Pickled Result.")

    # Save the results
    with open(pickle_path, 'wb') as outp:
        pickle.dump(output, outp, pickle.HIGHEST_PROTOCOL)

    return output, elapsed


# output.bg
# output.bc
# output.bunit

#bc: B on the centres of the cell faces.
#bg: B as a (weighted) averaged on grid points.
#bunit: Unit of the input map data.


def pixel_to_latlon(mag_map, header, fluxon_location):
    br = sunpy.map.Map(mag_map, header)
    shp = mag_map.shape
    n_lat = shp[0]
    lat_center = n_lat // 2
    f_lon = np.deg2rad(fluxon_location[:,0] *   br.meta['cdelt1'])
    f_lat = (fluxon_location[:,1]-lat_center) * br.meta['cdelt2']
    f_sgn = fluxon_location[:,2]
    n_flux = len(f_sgn)
    return f_lat, f_lon, f_sgn, n_flux


def get_fluxon_locations(floc_path, batch):
    fluxon_location = np.genfromtxt(floc_path)
    magnet, header = load_fits_magnetogram(batch=batch, ret_all=True)
    f_lat, f_lon, f_sgn, n_flux = pixel_to_latlon(magnet, header, fluxon_location)
    return f_lat, f_lon, f_sgn, n_flux

def plot_fluxon_locations(br_safe, cr, datdir, fits_path, reduce, force_plot=False, batch='fluxon', nwant=0, do_plot=False):
    # CL - DOUBLE CHECK THE LATITUTE DIRECTION HERE, and plot the polarity inversion line.
    floc_path = f"{datdir}/batches/{batch}/cr{cr}/floc/floc_cr{cr}_r{reduce}_f{nwant}.dat"
    f_lat, f_lon, f_sgn, n_flux = get_fluxon_locations(floc_path, batch)
    # fluxon_location = np.genfromtxt(floc_path)
    # magnet, header = load_fits_magnetogram(batch=batch, ret_all=True)
    # f_lat, f_lon, f_sgn, n_flux = pixel_to_latlon(magnet, header, fluxon_location)

    if not do_plot:
        return f_lat, f_lon, f_sgn

    the_dir = path.dirname(path.dirname(fits_path))
    fluxon_map_output_path = path.join(the_dir, f'r{reduce}_f{n_flux}_footpoints.pdf')

    top_dir = path.join(datdir,f"batches/{batch}/imgs/footpoints")
    fluxon_map_output_path_top = path.join(top_dir, f'cr{cr}_footpoints_r{reduce}_f{nwant}.pdf')

    fluxon_map_output_path_blank =      path.join(the_dir, "magnetograms", f'footpoints_cr{cr}_r{reduce}_blank.png')
    fluxon_map_output_path_blank_top =  path.join(top_dir, f'cr{cr}_footpoints_r{reduce}_blank.png')
    # fluxon_map_output_path2 = path.join(path.dirname(path.dirname(fits_path)), f'r{reduce}_f{n_flux}_footpoints_latlon.png')


    if not path.exists(top_dir):
        os.makedirs(top_dir)

    plot_paths = [
                fluxon_map_output_path_top,
                fluxon_map_output_path_blank_top,
                fluxon_map_output_path_blank, 
                ]

    print("\n\tPlotting Fluxon Footpoint Locations...", end="")

    need_plot = False
    for testpath in plot_paths:
        if not path.exists(testpath):
            need_plot = True
    if need_plot or force_plot:
        # br = sunpy.map.Map(datdir + 'hmi.Synoptic_Mr.polfil/hmi.synoptic_mr_polfil_720s.' + cr + '.Mr_polfil.fits')
        # fluxon_map_output_path = fits_path.replace('.fits', '_fluxons.png')
        # fluxon_map_output_path = path.join(path.dirname(path.dirname(fits_path)), path.basename(fits_path).replace('.fits', f'_{n_flux}_fluxons.png'))

        ## Print the Fluxon Map
        fig, ax = plt.subplots()

        magnet = br_safe.data
        # find the max and min of the magnetogram plot for use in setting the colormap, 
        sigma = 2
        mmean = np.nanmean(magnet)
        msig = np.nanstd(magnet)
        
        mvmin = mmean - sigma*msig
        mvmax = mmean + sigma*msig
        mvmin, mvmax = -500, 500

        # Plot the magnetogram
            # ax.imshow(br_safe.data, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), vmin=mvmin, vmax=mvmax,aspect='auto')


        # plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)

        # ax.imshow(br_safe.data, cmap='gray', interpolation=None, origin="lower", zorder=-10)
        # ratio = shp[1]/shp[0]
        # for (x, y, sig) in fluxon_location:
        #     color = 'red' if sig > 0 else "teal"
        #     ax.scatter(x, y, c=color, alpha=0.4)



        magimg = ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)
        shp=magnet.shape
        plt.axis('off')
        sz0=6 #inches
        ratio = shp[1]/shp[0]
        sz1=sz0*ratio #inches
        shp #pixels
        DPI = shp[1] / sz1 #pixels/inch
        fig.set_size_inches((sz1, sz0))
        plt.tight_layout()

        # plot a blank version of the map
        plt.savefig(fluxon_map_output_path_blank, bbox_inches='tight', dpi=4*DPI)

        # plot a blank version of the map
        magimg = ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)
        plt.savefig(fluxon_map_output_path_blank_top, bbox_inches='tight', dpi=4*DPI)

        # scatter the fluxons on top of the map
        magimg = ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower", aspect='auto', vmin=mvmin, vmax=mvmax)

        x, y, sig = zip(*fluxon_location)
        colors = ['red' if s > 0 else 'teal' for s in sig]
        ax.scatter(x, y, c=colors, alpha=0.4)

        # plot the fluxons scattered on top of the map
        # plt.savefig(fluxon_map_output_path, bbox_inches='tight', dpi=4*DPI)
        plt.savefig(fluxon_map_output_path_top, bbox_inches='tight', dpi=4*DPI)
        plt.close()

        ## Plot the fluxon map in lat/lon
    
        # fig, ax = plt.subplots()
        # ax.imshow(br_safe.data, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), vmin=mvmin, vmax=mvmax,aspect='auto')
        
        # # Plot the fluxons
        # colors = ['red' if s > 0 else 'teal' for s in f_sgn]
        # ax.scatter(f_lon, f_lat, c=colors, alpha=0.4)

        # # for (long, lat, sig) in zip(f_lon, f_lat, f_sgn):
        # #     color = 'red' if sig > 0 else "teal"
        # #     ax.scatter(long, lat, color=color)
        # # # rewrite the previous for loop without a loop


        # # fluxon_map_output_path2 = path.join(path.dirname(path.dirname(fits_path)), path.basename(fits_path).replace('.fits', f'_latlon_{n_flux}_fluxons.png'))
        # plt.axis('off')
        # fig.set_size_inches((sz1, sz0))
        # plt.tight_layout()
        # plt.savefig(fluxon_map_output_path2, bbox_inches='tight', dpi=4*DPI)
        # plt.close()
        # plt.show()
        print("Success! Files saved to: ", end="\n\t\t")
    else:
        print("Skipped! Files already exist:", end="\n\t\t")
    for testpath in plot_paths:
        print(shorten_path(testpath, 5), end="\n\t\t")
    return f_lat, f_lon, f_sgn

from tqdm import tqdm
import timeout_decorator
from pfsspy import tracing
import astropy.constants as const




def trace_lines(output, f_lon, f_lat, f_sgn, open_path, closed_path):


    # Note that this code was originally developed using an older version of the pfsspy code - improvements look to be able to handle bundles of fieldlines, which would greatly simplify some of the code. Do this.



    # Capture an array of starting and endpoints for open fieldlines 2x(r, lat, lon) + sgn of starting point
    fl_open = np.zeros([1,5])
    flnum_open = 0

    # Capture an array of starting and endpoints for closed fieldlines 2x(r, lat, lon) + sgn of starting point
    fl_closed = np.zeros([1,5])
    flnum_closed = 0

    r0 = 1.01 * const.R_sun



    skip_num = 0
    timeout_num = 0
    tracer = tracing.PythonTracer()
    for i, coords in enumerate(tqdm(zip(f_lon, f_lat), desc="Tracing Field Lines", total=len(f_lat))):
        try:
            output, fl_open, fl_closed, flnum_open, flnum_closed, tracer = trace_each(coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer, r0, f_sgn)
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
        print(f"\n\n\t\tSome iterations failed. Timed-out: {timeout_num} ({t_perc:0.2f}%), ValueError: {skip_num} ({s_perc:0.2f}%)\n")


    print(f"\n\tOpen Lines: {flnum_open+1}, Closed Lines: {flnum_closed}, Failures: {skip_num+timeout_num}, Total Good: {flnum_open+flnum_closed}")
    # print(f"\n\tOpen Lines: {len(fl_open)}, Closed Lines: {len(fl_closed)}, Failures: {skip_num+timeout_num}")
        
    fl_open = fl_open[1:]
    fl_closed = fl_closed[1:]


    # Output is flnum, polarity, latitude, longitude, radius
    # Save these coordinates to file
    print("\n Saving Fluxons...", end="")
    np.savetxt(open_path, fl_open)
    np.savetxt(closed_path, fl_closed)
    print("Success!")

    return fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed



@timeout_decorator.timeout(5)
def trace_each(coords, i, output, fl_open, fl_closed, flnum_open, flnum_closed, tracer, r0, f_sgn):
    (this_flon, this_flat) = coords
    coord_frame = output.coordinate_frame
    # THE ARCSIN IS A TEST
    x0 = SkyCoord(this_flon * u.rad, np.arcsin(this_flat) * u.rad, r0, frame=coord_frame)
    # x0 = SkyCoord(this_flon * u.rad, this_flat * u.rad, r0, frame=coord_frame)
    # x0_lat = x0.lat
    # x0_lon = x0.lon
    # xx_lat = np.rad2deg(np.arcsin(this_flat))
    # xx_lat2 = np.rad2deg(this_flat)
    # xx_lon = np.rad2deg(this_flon)
    fl = output.trace(tracer, x0)
    fl = fl.field_lines[0]
    if fl.is_open:
        fl_pol = fl.polarity
        fl = fl.coords.spherical
        fl_lats = fl.lat.value
        fl_lons = fl.lon.value
        fl_rads = (fl.distance/const.R_sun).value

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter(fl_lats, fl_lons, fl_rads)
        # plt.show()

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