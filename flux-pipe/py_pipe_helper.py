import sunpy
import sunpy.io
import sunpy.coordinates
# import sunpy.net
from sunpy.net import Fido, attrs as a
# Fido = sunpy.net.Fido
import drms
import os
import glob
import sunpy.coordinates.frames as frames
import astropy.units as u
default_email = "chris.gilly@colorado.edu"
from pathlib import PosixPath
# import subprocess
from astropy.nddata import block_reduce
from astropy.io import fits
import numpy as np
import os
import os.path
import sys
# import ADAPTClient
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt # Import libraries
# from sunpy.net.dataretriever import GenericClient


def add_paths(flux_pipe_dir):
    # Path to the PDL script
    pdl_script_path = flux_pipe_dir + "magnetogram2wind.pdl"
    os.chdir(flux_pipe_dir)
    # Get the plotscript directory path
    plot_dir = os.path.abspath(os.path.join(flux_pipe_dir, "plotting"))
    sys.path.append(plot_dir)
    return pdl_script_path

# Magnetogram things

def make_mag_dir(datdir, ADAPT=False):
    mag_dir = os.path.join(datdir, "magnetograms")

    if ADAPT:
        mag_dir = os.path.join(mag_dir, "ADAPT")

    if not os.path.exists(mag_dir):
        os.makedirs(mag_dir)
    return mag_dir

def get_magnetogram_file(cr=None, date=None, datdir=None, email=None, force_download=False, reduce = False):
    """
    Function to grab HMI data.

    Args:
        cr (int): Carrington rotation number.
        date (str): Date in YYYY-MM-DD format.
        data_dir (str): Optional directory where data will be stored. If not specified,
            default directories will be used.

    Returns:
        None
    """

    # Set the download account
    try:
        jsoc_email = email or os.environ["JSOC_EMAIL"]
    except KeyError:
        jsoc_email = default_email


    # Set the Carrington rotation number
    if cr is not None:
        CR = cr
        date = sunpy.coordinates.sun.carrington_rotation_time(CR)
    elif date is not None:
        CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))
    else:
        raise ValueError("Must specify either cr or date!")



    mag_dir = make_mag_dir(datdir)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"(py) Getting Magnetogram for CR{CR}, from {date}...")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    import pathlib
    hmi_object = pathlib.Path(mag_dir)
    file_list = list(hmi_object.iterdir())
    print("\tSearching for file...")
    found_file = False
    for file in file_list:
        if str(CR)+"_r1_" in str(file):
            print(f"\t\tFound '{os.path.basename(file)}' in '{shorten_path(mag_dir)}'")
            found_file = True
            break

    if found_file:
        if force_download:
            print("\tForcing redownload...")
        else:
            small_path = reduce_mag_file(file, reduce, force=force_download)
            return file, small_path


    c = drms.Client()
    # Generate a search
    crot = a.jsoc.PrimeKey('CAR_ROT', str(CR))
    res = Fido.search(a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), crot,
                    a.jsoc.Notify(jsoc_email))

    # Once the query is made and trimmed down...
    big_path = os.path.join(mag_dir, f"CR{CR}_r1_hmi.fits")
    # hmi_path = hmidat+'/{file}.fits'


    print("\tDownloading HMI from JSOC...")
    out = Fido.fetch(res, path=mag_dir)
    hmi_path_out = out[0]
    os.rename(hmi_path_out, big_path)
    print(f"\n\tSaved to {big_path}\n")


    small_path = reduce_mag_file(big_path, reduce, force=force_download)
    return big_path, small_path


def get_ADAPT_file(cr=None, date=None, datdir=None, email=None, force_download=False, reduce = False, method=2):
    """
    Function to grab ADAPT data.

    Args:
        cr (int): Carrington rotation number.
        date (str): Date in YYYY-MM-DD format.
        data_dir (str): Optional directory where data will be stored. If not specified,
            default directories will be used.

    Returns:
        None
    """

    print(reduce)
    print("_____________________")

    ## Parse the Dates
    # Set the Carrington rotation number
    if cr is not None:
        CR = cr
        date = sunpy.coordinates.sun.carrington_rotation_time(CR)
        # date_end = sunpy.coordinates.sun.carrington_rotation_time(CR+1)
    elif date is not None:
        CR = int(sunpy.coordinates.sun.carrington_rotation_number(date))
        # date_end = sunpy.coordinates.sun.carrington_rotation_time(CR+1)
    else:
        raise ValueError("Must specify either cr or date!")

    date_end = date + (1.9999999 * u.hour)

    # Format the Display Dates
    tstring_display =r"%H:%M:%S %m-%d-%Y"
    display_date    = date.strftime(tstring_display)
    display_date_end = date_end.strftime(tstring_display)

    # Format the Search Dates
    tstring = r"%Y-%m-%dT%H:%M:%S"
    get_date    =date.strftime(tstring)
    get_date_end=date_end.strftime(tstring)

    # Make the directory
    mag_dir = make_mag_dir(datdir, ADAPT=True)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"(py) Getting ADAPT Magnetogram(s) for CR{CR}, from {display_date} to {display_date_end}...")
    print(  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    print("\tChecking for file...")
    found_file = False

    date_obs = [date, date + (1.9999999 * u.hour), date - (1.9999999 * u.hour)]

    dates = [xx.strftime("%Y%m%d") for xx in date_obs]

    import pathlib
    path_obj = pathlib.Path(mag_dir)
    file_list = list(path_obj.iterdir())

    for file in file_list:
        look = str(CR)+"_r1_"
        file_string = str(file)
        if look in file_string:
            found_file = True
            break
            for ii, dt in enumerate(dates):
                if dt in str(file):
                    print(f"\t\tFound '{os.path.basename(file)}' in '{shorten_path(mag_dir)}'")
                    print("\t\tDate: ", date_obs[ii].strftime(tstring_display))
                    found_file = True
                    break
            if found_file:
                break

    if found_file:
        if force_download:
            print("\tForcing redownload...")
        else:
            print("\tSkipping Download!\n")
            # small_path = None # reduce_mag_file(file, reduce, force=force_download)
            # mean_path = format_ADAPT_file(file, reduce, force=force_download)
            the_path = format_ADAPT_file(file, method=method, force=force_download)

            return file, the_path
    else:
        print("\t\tNo file found!")

    print("\n\tSearching FIDO for ADAPT Map...\n")
    from ADAPTClient import ADAPTLngType
    LngType = '0' # 0 is carrington, 1 is central meridian
    print("VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV")
    res = Fido.search(a.Instrument('adapt'), a.Time(get_date, get_date_end), ADAPTLngType(LngType))
    print(res)
    print("\tDownloading ADAPT map...\n")
    out = Fido.fetch(res, path=mag_dir)
    assert len(out) == 1, "More than one file found!"
    if len(out) ==1: print("\tSuccess!")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")


    file_all = str(out[0])
    # file_end = file_all.split('adapt')[1].split('.')
    # file_end = file_all.split('adapt')[1].split('.')[1:]
    big_path = os.path.join(mag_dir, f"CR{CR}_r1_adapt.fts.gz")


    os.rename(file_all, big_path)
    print(f"\n\t\tSaved to {big_path}")

    # mean_path = format_ADAPT_file(big_path, reduce, force=force_download)
    the_path = format_ADAPT_file(big_path, method=method, force=force_download)

    return big_path, the_path

def expose_adapt_fits(adapt_fits):
    return  [print(x, " : ", adapt_fits[0].header[x]) for x in adapt_fits[0].header if not 'KEYCOMMENTS' in x]

def print_adapt_maps(adapt_maps):
    for map in adapt_maps:
        print(map)

def plot_all_adapt(adapt_maps):
    from matplotlib import gridspec
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(4, 3, figure=fig)
    for i, a_map in enumerate(adapt_maps):
        ax = fig.add_subplot(gs[i], projection=a_map)
        a_map.plot(axes=ax, cmap='bwr', vmin=-2, vmax=2,
                title=f"Realization {1+i:02d}")

    # adapt_maps.plot()
    plt.tight_layout(pad=5, h_pad=2)
    plt.show(block=True)

def plot_mean_adapt(mean_adapt):
    map_mean = np.nanmean(mean_adapt)
    map_std = np.std(mean_adapt)
    plt.imshow(mean_adapt, vmin = map_mean - 2*map_std, vmax=map_mean + 2*map_std)
    plt.show(block=True)

def format_ADAPT_file(filename, method='mean', force=False):
    import sunpy.io
    import sunpy.map

    if method == 'mean':
        out_file_name = str(filename).replace("_r1_", "_rmean_")
    else:
        out_file_name = str(filename).replace("_r1_", f"_rf{method}_")

    out_file_name = out_file_name.replace(".fts.gz", ".fits")
    # print(out_file_name)

    if os.path.exists(out_file_name) and not force:
        print("\tFile already formatted!")
        print("\t\t", shorten_path(out_file_name), "\n")
        print("\t\t\t```````````````````````````````\n\n")

        return out_file_name

    print("\n\tFormatting ADAPT file...", end='')
    adapt_fits = sunpy.io.read_file(filename)
    main_header = adapt_fits[0].header
    main_data = adapt_fits[0].data
    main_header['DATE-AVG'] = main_header['MAPTIME']

    with fits.open(filename, ignore_missing_simple=True) as hdul:
        hdul.verify('silentfix')
        header2 = hdul[0].header



    if method == 'mean':
        data_header_pairs = [(map_slice, main_header) for map_slice in main_data]
        adapt_maps = sunpy.map.Map(data_header_pairs, sequence=True)
        adapt_cube = np.asarray([the_map.data for the_map in adapt_maps])
        output_map = np.nanmean(adapt_cube, axis=0)

        if False:
            # Lots of Plots
            plot_all_adapt(adapt_maps)
            expose_adapt_fits(adapt_fits)
            print_adapt_maps(adapt_maps)
            plot_mean_adapt(mean_adapt)

    elif type(method) == int:
        adapt_map = sunpy.map.Map((main_data[method], main_header))
        output_map = np.asarray(adapt_map.data)
    else:
        assert False, "Method not recognized!"

    useheader = fix_header_ADAPT(header2, output_map)
    fits.writeto(out_file_name, output_map, useheader, overwrite=True)

    print("Success!")
    print("\t\tSaved to", out_file_name, "\n")


    return out_file_name


def reduce_mag_file(mag_file, reduction=3, force=False):
    """Reduces the size of a magnetogram FITS file by a given factor."""
    small_file = PosixPath(str(mag_file).replace("_r1_", f"_r{reduction}_"))
    # reduce the FITS image
    print(f"\tReducing image size by a factor of {reduction}...", end="")
    if not os.path.exists(small_file) or force:
        small_file = reduce_fits_image(mag_file, small_file, target_resolution=None, reduction_amount=reduction)
        # print("Success!\n")
    else:
        print("Skipped! Reduced file already exists:")
        # print("\t\t", shorten_path(str(small_file), 2))
        ### WORKING HERE
        print(f"\t\tFound '{os.path.basename(small_file)}' in '{shorten_path(os.path.dirname(small_file))}'")
        print("\n\t\t\t```````````````````````````````\n \n\n")

    return small_file

def read_fits_data(fname):
    """Reads FITS data and fixes/ignores any non-standard FITS keywords."""
    hdulist = fits.open(fname, ignore_missing_simple=True)
    hdulist.verify('silentfix+warn')
    return hdulist

def reduce_fits_image(fits_path, small_file, target_resolution=None, reduction_amount=None, func=np.nansum):
    """
    Open a FITS file, reduce the size of the image using astropy's block_reduce
    function, and save a new copy of the FITS file with the smaller image in the
    same directory as the original.

    :param fits_path: str, path to the FITS file
    :param target_resolution: int, optional, target resolution in pixels, if
                              specified reduction_amount is ignored
    :param reduction_amount: int, optional, amount to reduce the size of the
                              image, if target_resolution is not specified
    :param func: numpy function, optional, function to use for the reduction
                              defaults to np.nanmean
    """
    print(f"\n\tReducing {fits_path}...")
    # Open the FITS file and read the data
    with fits.open(fits_path, ignore_missing_simple=True) as hdul:
        hdul.verify('silentfix')
        data = hdul[0].data
        if data is None:
            data = hdul[1].data

        current_resolution = max(data.shape)
        print("\tOriginal Shape: ", data.shape)

        # Calculate the reduction amount if target resolution is specified
        if target_resolution is not None:
            reduction_amount = int(np.ceil(current_resolution / target_resolution))

        # Raise an error if neither target_resolution nor reduction_amount is specified
        elif reduction_amount is None:
            raise ValueError("Either target_resolution or reduction_amount must be specified.")


        before_sum = np.sum(data)
        small_image = block_reduce(data, reduction_amount, func)
        after_sum = np.sum(small_image)
        if not np.isclose(before_sum, after_sum):
            print("\tREDUCTION WARNING: \n\tSum before:    ", before_sum, "\n\tSum after:     ", after_sum)

        try:
            hdul[0].header["DATE"]
            useheader = hdul[0].header
        except KeyError:
            useheader = hdul[1].header

        del useheader['BLANK']
        useheader = fix_header(useheader, small_image)
        # small_file = fits_path.replace('_r1_', f'_r{reduction_amount}_')

        # del useheader['BLANK']
        # useheader['DATAMIN'] = np.min(small_image)
        # useheader['DATAMAX'] = np.max(small_image)
        # useheader['BZERO'] = 0
        # useheader['BSCALE'] = 1

        # useheader['CDELT1'] = 360 / small_image.shape[1]  ## DEGREES
        # useheader['CDELT2'] = np.deg2rad(360 / (small_image.shape[0] * np.pi)) #RADIANS

        print("\tFinal Shape:    ", small_image.shape)

        print("\tSaving  ", small_file)
        fits.writeto(small_file, small_image, useheader, overwrite=True)

        # plot_raw_magnetogram(fits_path, data, small_image)

        print("    Reduction Complete!\n")
    print("```````````````````````````\n")

    return small_file

def fix_header(useheader, image):

    useheader['DATAMIN'] = np.min(image)
    useheader['DATAMAX'] = np.max(image)
    useheader['BZERO'] = 0
    useheader['BSCALE'] = 1

    useheader['CDELT1'] = 360 / image.shape[1]  ## DEGREES
    useheader['CDELT2'] = np.deg2rad(360 / (image.shape[0] * np.pi)) #RADIANS
    return useheader


def fix_header_ADAPT(useheader, image):
    useheader['DATAMIN'] = np.min(image)
    useheader['DATAMAX'] = np.max(image)
    useheader['BZERO'] = 0
    useheader['BSCALE'] = 1

    # import pdb; pdb.set_trace()
    useheader['CDELT1'] = 360 / image.shape[1]  ## DEGREES
    useheader['CDELT2'] = 360 / (image.shape[0] * np.pi) ## DEGREES
    return useheader

def plot_raw_magnetogram(fits_path, data, small_image):
    # Save the high resolution image as a grayscale PNG
    plt.axis('off')
    high_res_output_path = fits_path.replace('.fits', '.png')
    fig = plt.gcf()
    shp = data.shape
    dmin = np.nanmin(data)
    dmax = np.nanmax(data)
    dmean = np.nanmean(data)
    dsig = np.nanstd(data)
    thresh = 3
    vmin = dmean - thresh*dsig
    vmax = dmean + thresh*dsig
    plt.imshow(data, cmap='gray', vmin=vmin, vmax=vmax)

    ratio = shp[1]/shp[0]
    sz0=6 #inches
    sz1=sz0*ratio #inches
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.savefig(high_res_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()

    # Save the low resolution image as a grayscale PNG
    plt.imshow(small_image, cmap='gray', vmin=vmin, vmax=vmax)
    plt.axis('off')
    low_res_output_path = fits_path.replace('.fits', '_small.png')
    fig = plt.gcf()
    shp = small_image.shape
    ratio = shp[1]/shp[0]
    sz0=6 #inches
    sz1=sz0*ratio #inches
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.savefig(low_res_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()

def load_fits_magnetogram(datdir = "/Users/cgilbert/vscode/fluxon-data/", batch="fluxon", bo=2, bn=2, ret_all=False):
    """Loads a magnetogram from a FITS file."""
    fname = load_magnetogram_params(datdir)[2].replace("/fluxon/", f"/{batch}/").replace(f"_{bo}_", f"_{bn}_")
    fits_path = datdir + fname
    try:
        hdulist = read_fits_data(fits_path)
    except FileNotFoundError as e:
        hdulist = read_fits_data(fname)
    brdat = hdulist[0].data
    header= hdulist[0].header
    brdat = brdat - np.mean(brdat)
    if ret_all:
        return brdat, header
    else:
        return brdat


# File I/O and pathing
def write_magnetogram_params(datdir, cr, file_path, reduction, adapt=False):
    """Writes the magnetic_target.params file for a given CR and reduction amount."""
    # write the parameter file
    adapt_str = 1 if adapt else 0
    reduce = 'A' if adapt else reduction
    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'w') as fp:
        fp.write("## CR_int, Filename_str, Adapt_bool, Doplot_bool, reduction ##\n")
        fp.write(str(cr)+"\n")
        fp.write(str(file_path)+"\n")
        fp.write(str(adapt_str)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(reduce))

def load_magnetogram_params(datdir):
    """Reads the magnetic_target.params file and returns the parameters."""
    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'r') as fp:
        hdr = fp.readline().rstrip()
        cr = fp.readline().rstrip()
        fname = fp.readline().rstrip()
        adapt = int(fp.readline().rstrip())
        doplot = int(fp.readline().rstrip())
        reduce = fp.readline().rstrip()
    return (hdr, cr, fname, adapt, doplot, reduce)

def find_file_with_string(directory, search_string):
    """Searches a directory for a file containing a given string."""
    for file_name in os.listdir(directory):
        if search_string in file_name:
            return os.path.join(directory, file_name)
    return None

def shorten_path_old(start_path, levels=1):
    start_parts = start_path.split('/')
    out_parts = start_parts[-levels:] if levels > 0 else start_parts
    out_string = '/'.join(out_parts) if out_parts else start_path
    return "DATAPATH/" + out_string

def shorten_path(string, __=None):
    datapath = os.getenv("DATAPATH")
    if datapath:
        return string.replace(datapath, "$DATAPATH")
    else:
        return string



def get_fixed_coords(phi0, theta0):
    ph0, th0 = phi0+np.pi, np.sin(-(theta0-(np.pi/2)))
    return ph0, th0



# def find_hilbert_footpoints(batchdir, cr, want_points=1000, reduction=3, force=False):
#     """Finds the hilbert footpoints for a given CR and reduction amount."""
#     N_fluxons = want_points

#     # get the data directory
#     flocdir = os.path.join(batchdir, f"cr{cr}", "flocs")
#     flocpath = os.path.join(flocdir, f'foot_locs_cr{cr}_{want_points}.dat')
#     mag_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(flocdir))))
#     mag_path = os.path.join(mag_dir, f"magnetograms/CR{cr}_r{reduction}_hmi.fits")
#     out_path = None

#     if not os.path.exists(flocdir):
#         os.makedirs(flocdir)

#     if not os.path.exists(flocpath) or force:

#         # with fits.open(mag_path, ignore_missing_simple=True) as hdul:
#         #     print("\n\t**Reading Magnetogram...")
#         #     hdul.verify('silentfix')
#         #     data = hdul[0].data
#         #     if data is None:
#         #         data = hdul[1].data

#     #     # if adapt:
#     #     #     smag = smag[:,:,2].squeeze()

#             # print("\n")
#         print("\n\t**Placing {0} footpoints using Hilbert Curves...\n".format(N_fluxons))

#         tolerance = N_fluxons * 0.05
#         step = N_fluxons * 0.1
#         width = 100
#         most = N_fluxons + width
#         least = N_fluxons - width
#         iterate = 1
#         try_fluxons = int(N_fluxons * 1.25)
#         count = 0
#         maxreps = 30
#         codepath = os.path.join(os.path.dirname(mag_dir), "fluxon-mhd/pdl/PDL/")
#         hilbertpath = os.path.join(codepath, "fluxon_placement_hilbert.pdl")
#         # os.system(f"cd {codepath}")

#         while iterate and count < maxreps:
#             try_fluxons = int(try_fluxons)
#             result = subprocess.run(["perl", hilbertpath, mag_path, str(try_fluxons), flocpath])
#             # floc = fluxon_placement_hilbert(smag, try_fluxons)

#             exit()
#             N_actual = len(floc) / 3
#             count += 1

#             print("\t  Placing footpoints:: iter: {0}/{1}. Wanted: {2}, Tried: {3}, Placed: {4}...".format(count, maxreps, N_fluxons, try_fluxons, N_actual))

#             distance = N_actual - N_fluxons
#             if abs(distance) > tolerance:
#                 factor = 1 + abs(distance) / (N_fluxons * 2)
#                 step *= factor

#                 if distance > 0:
#                     try_fluxons -= step
#                 else:
#                     try_fluxons += step

#                 iterate = 1
#                 print("Retrying!\n")
#             else:
#                 print("Success!\n")
#                 iterate = 0

#     #     print("\n")
#     #     print("\t**Writing Result to Disk...")
#     #     np.savetxt(out_path, floc.transpose())
#     # else:
#     #     floc = np.loadtxt(out_path)
#     #     floc = floc.transpose()
#     #     N_actual = int(len(floc))
#     #     print("Skipped! Found {0} footpoints on Disk.".format(N_actual))

#     print("\n")



#     # # get the magnetogram file name
#     # mag_file = find_highest_numbered_file(datdir, "br_r1")
#     # # reduce the FITS image
#     # small_file = reduce_mag_file(mag_file, reduction=reduction, force=force)
#     # # write the parameter file
#     # write_magnetogram_params(datdir, cr, small_file, reduction)
#     # # run the hilbert footpoint finder
#     # run_hilbert_footpoint_finder(datdir, cr, reduction=reduction, force=force)
#     # # return the hilbert footpoint file
#     # return find_highest_numbered_file(datdir, "floc_cr"+str(cr))

#     # import os
#     # import numpy as np

#     # print("\n\n**Processing Magnetograms to get Footpoints...")

#     # no_floc = 0
#     # out_dir = os.path.join(datdir, batch_name, "cr" + str(cr), "floc")



# if __name__ == "__main__":
#     get_magnetogram_file(force_download=False, datdir="/Users/cgilbert/vscode/fluxon-data/")








  # print("\nConfiguring Download...")


    # # Specify any directories
    # if datdir is not None:
    #     if cr is not None:
    #         datdir += f"/{batch}/cr{cr}/mag"
    #     hmidat = os.path.expanduser(os.path.join(datdir, 'hmi.Synoptic_Mr.polfil'))
    #     mdidat = os.path.expanduser(os.path.join(datdir, 'mdi.Synoptic_Mr.polfil'))
    # else:
    #     hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil')
    #     mdidat = os.path.expanduser('~/data/mdi.Synoptic_Mr.polfil')

    # HMI data

    # # Sort out the last downloaded rotation
    # crfiles = glob.glob(hmidat+'*.fits')
    # crfiles.sort()
    # crlist = [int(i[-19:-15]) for i in crfiles]

    # Specify requested rotations

    # if cr

    # if cr is None and date is None:
    #     if len(crlist) == 0:
    #         cr0 = 2096
    #     else:
    #         cr0 = max(crlist) + 1
    #     cr1 = int(sunpy.coordinates.sun.carrington_rotation_number(t='now'))

    #     if (cr0 - 1) == cr1:
    #         print('what?')
    # else:
    #     if cr is not None:
    #         cr0 = cr
    #     else:
    #         cr0 = int(sunpy.coordinates.sun.carrington_rotation_number(date))

    #     cr1 = cr0


        # MDI data

    # # Grab MDI
    # if 2104 >= cr >= 1911:
    #     print("\nDownloading MDI from JSOC...")
    #     os.system('mkdir ' + mdidat)
    #     os.chdir(mdidat)
    #     mdi_file = "synop_Mr_0.polfil.{}.fits".format(cr1)
    #     mdi_path = os.path.join(mdidat,mdi_file)
    #     address = f"http://soi.stanford.edu/magnetic/synoptic/carrot/M_Corr/{mdi_file}"
    #     command = 'curl -O "{}"'.format(address)
    #     os.system(command)
    #     # os.system('curl -O "http://soi.stanford.edu/magnetic/synoptic/carrot/M_Corr/synop_Mr_0.polfil.[1911-2104].fits"')
    # else:
    #     mdi_path = None
    #     print("\n !! No MDI data available for this time period !!\n")

    # print("\tDownload Complete!\n")



    # ADAPTFileType
    # ADAPTLngType
    # ADAPTInputSource
    # ADAPTDataAssimilation
    # ADAPTResolution
    # ADAPTVersionYear
    # ADAPTVersionMonth
    # ADAPTEvolutionMode
    # ADAPTHelioData
    # ADAPTMagData



    # print(baseurl)
    # from sunpy.net import Scraper
    # scraper = Scraper(baseurl, regex=True)
    # import pdb; pdb.set_trace()

    # scraper = Scraper(baseurl, pattern, regex=True)

    # @classmethod
    # def register_time(cls):
    #     map_time = a.TimeParameter('map_time', a.Time)
    #     # print(map_time)
    #     return map_time
        # def _get_url_for_timerange(self, timerange, **kwargs):
    #     """
    #     Returns list of URLs corresponding to value of input timerange.
    #     """
    #     from sunpy.net import scraper
    #     # import pdb; pdb.set_trace()
    #     time_scraper = scraper.Scraper(self.baseurl + self.pattern)
    #     thing = time_scraper.filelist(timerange)
    #     # print(thing)
    #     1/0
    #     return thing



    # @classmethod
    # def register_values(cls):
    #     adict = {attrs.Instrument: [('ADAPT', 'ADvanced Adaptive Prediction Technique.')],
    #              attrs.Source: [('NSO', 'National Solar Observatory.')],
    #              attrs.Provider: [('GONG', 'Global Oscillation Network Group.')],
    #              attrs.Time: []
    #              }
    #     return adict

# Fido.register_client(ADAPTClient)







    # print("ASDF")
    # 1/0

    # Set the download account


        # baseurl = url + r'%Y'
  # pattern = r'{year:4d}/adapt{file_type:d}{lng_type:d}{input_source:d}{data_assimilation:d}{resolution:d}_{version_yr:d}{version_month:l}{realizations:d}_{map_time}_{evolution_mode:l}{days_since_last_obs:d}{hours_since_last_obs:d}{minutes_since_last_obs:d}{seconds_since_last_obs:d}{helio_data:l}{mag_data:d}.fts.gz'
    # pattern = r'{file_type:s}{lng_type:s}{input_source:s}{data_assimilation:s}{resolution:s}_{version_yr:s}{version_month:s}{realizations:d}_{map_time}_{evolution_mode:s}{days_since_last_obs:d}{hours_since_last_obs:d}{minutes_since_last_obs:d}{seconds_since_last_obs:d}{helio_data:s}{mag_data:s}.fts.gz'
    # pattern = "{}"
    # (\d){5}_(\d){12}_(\w){1}(\d){8}(\w){1}(\d){1}\.fts\.gz
    # pattern = r'{}adapt{ADAPTFileType:1d}{ADAPTLngType:1d}{ADAPTInputSource:1d}{ADAPTDataAssimilation:1d}{ADAPTResolution:1d}_{ADAPTVersionYear:2d}{ADAPTVersionMonth:1s}{realizations:03d}_{map_time:12d}_{ADAPTEvolutionMode:1s}{days_since_last_obs:2d}{hours_since_last_obs:2d}{minutes_since_last_obs:2d}{ADAPTHelioData:1s}{ADAPTMagData:1d}.fts.gz'
    # pattern = '{}adapt{ADAPTFileType:1d}{ADAPTLngType:1d}{ADAPTInputSource:1d}{ADAPTDataAssimilation:1d}{ADAPTResolution:1d}_{ADAPTVersionYear:2d}{ADAPTVersionMonth:1l}{realizations:3d}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}_{ADAPTEvolutionMode:1l}{days_since_last_obs:2d}{hours_since_last_obs:2d}{minutes_since_last_obs:2d}{ADAPTHelioData:1l}{ADAPTMagData:1d}.fts.gz'
    # pattern = '{}adapt{ADAPTFileType:1d}{ADAPTLngType:1d}{ADAPTInputSource:1d}{ADAPTDataAssimilation:1d}{ADAPTResolution:1d}_{ADAPTVersionYear:2d}{ADAPTVersionMonth:1l}{realizations:3d}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}_{ADAPTEvolutionMode:1l}{days_since_last_obs:2d}{hours_since_last_obs:2d}{minutes_since_last_obs:2d}{ADAPTHelioData:1l}{ADAPTMagData:1d}.fts.gz'
    # pattern = 'adapt{:5d}_{:2d}{:1s}{:3d}_{:12d}_{:1s}{:8d}{:1s}{:1d}.fts.gz'
    # pattern = '{}adapt{:5d}_{:2d}{:1l}{:3d}_{:12d}_{:1l}{:8d}{:1l}{:1d}.fts.gz'


# from sunpy.net.dataretriever.client import GenericClient
# class ADAPTClient(GenericClient):

#     baseurl = ("https://gong.nso.edu/adapt/maps/gong/"
#                )

#     baseurl = (r'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/'
#                 r'L0CS/SpWx/%Y/%Y%m%d_EVE_L0CS_DIODES_1m.txt')
#     pattern = '{}/SpWx/{:4d}/{year:4d}{month:2d}{day:2d}_EVE_L{Level:1d}{}'

#     @classmethod
#     def register_values(cls):
#         from sunpy.net import attrs
#         adict = {attrs.Instrument: [('EVE', 'Extreme ultraviolet Variability Experiment, which is part of the NASA Solar Dynamics Observatory mission.')],
#                 attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
#                 attrs.Source: [('SDO', 'The Solar Dynamics Observatory.')],
#                 attrs.Provider: [('LASP', 'The Laboratory for Atmospheric and Space Physics.')],
#                 attrs.Level: [('0', 'EVE: The specific EVE client can only return Level 0C data. Any other number will use the VSO Client.')]}
#         return adict


# import re


    # baseurl = (r'https://gong.nso.edu/adapt/maps/gong/'
    #            r'adapt[ZXABR_CCEFFF]_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}_[TIIJJKKLLGQ].fts')
    # pattern = r'adapt(?P<modeltag>[ZXABR_CCEFFF])_(?P<datetime>{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d})_(?P<lastobstag>[TIIJJKKLLGQ]).fts'

    # baseurl = r'https://gong.nso.edu/adapt/maps/gong/' '\d{4}/adapt\d{5}_\d{2}[a-z]\d{3}_\d{12}_[a-z]\d{8}[a-z]\d{1}.fts.gz'
    # pattern = r'{}/{:d}/adapt{:d}_{:d}{:l}{:d}_{:d}{:d}{:d}{:d}_{:l}{:d}{:d}{:d}{:d}{:l}{:d}.fts.gz'
    # https://gong.nso.edu/adapt/maps/gong/2010/adapt40311_03k012_201005010200_i00100900n1.fts.gz

    # baseurl = r'https://gong.nso.edu/adapt/maps/gong/\d{4}/adapt4(\d{1})(\d{1})(\d{1})1_(\w{2})(\w{1})(\d{3})_\d{4}\d{2}\d{2}\d{2}\d{2}_i(\d{2})(\d{2})(\d{2})(\d{2})(\w{1})(\d{1}).fts.gz'
    # pattern = (r'https://gong.nso.edu/adapt/maps/gong/%Y/adapt4{LNGTYPE:1d}{InputSource:1d}{DataAssimilation:1d}1_{CodeVersionYear:2s}{CodeVersionMonth:1s}{Realization:3s}_%Y%m%d%H%M_i{DaysSinceLastObs:2d}{HoursSinceLastObs:2d}{MinutesSinceLastObs:2d}{SecondsSinceLastObs:2d}{HelioseismicData:1s}{MagData:1d}.fts.gz')

    # pattern = r'{year:4d}/adapt{file_type:d}{lng_type:d}{input_source:d}{data_assimilation:d}{resolution:d}_{version_yr:d}{version_month:l}{realizations:d}_{map_time_yr:d}{map_time_month:d}{map_time_day:d}{map_time_hour:d}{map_time_min:d}_{evolution_mode:l}{days_since_last_obs:d}{hours_since_last_obs:d}{minutes_since_last_obs:d}{seconds_since_last_obs:d}{helio_data:l}{mag_data:d}.fts.gz'
    # print(baseurl)
    # print(pattern)


        # from pfsspy.sample_data import get_adapt_map
    # adapt_fname = get_adapt_map()