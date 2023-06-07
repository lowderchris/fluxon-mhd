# Import libraries
import sunpy
import sunpy.io
import sunpy.coordinates
from sunpy.net import Fido, attrs as a
import drms
import os
import glob
import sunpy.coordinates.frames as frames
import astropy.units as u
default_email = "chris.gilly@colorado.edu"
from pathlib import PosixPath
import subprocess

def make_mag_dir(datdir):
    mag_dir = os.path.join(datdir, "magnetograms")
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
    print(f"Getting Magnetogram for CR{CR}, from {date}...")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    import pathlib
    hmi_object = pathlib.Path(mag_dir)
    file_list = list(hmi_object.iterdir())

    found_file = False
    for file in file_list:
        if str(CR)+"_r1_" in str(file):
            print(f"\tFound '{os.path.basename(file)}' in '{mag_dir}'")
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


from astropy.nddata import block_reduce
from astropy.io import fits
import numpy as np


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
        
        # small_file = fits_path.replace('_r1_', f'_r{reduction_amount}_')
        try:
            hdul[0].header["DATE"]
            useheader = hdul[0].header
        except KeyError:
            useheader = hdul[1].header

        del useheader['BLANK']
        useheader['DATAMIN'] = np.min(small_image)
        useheader['DATAMAX'] = np.max(small_image)
        useheader['BZERO'] = 0
        useheader['BSCALE'] = 1

        useheader['CDELT1'] = 360 / small_image.shape[1]  ## DEGREES
        useheader['CDELT2'] = np.deg2rad(360 / (small_image.shape[0] * np.pi)) #RADIANS

        print("\tFinal Shape:    ", small_image.shape)

        print("\tSaving  ", small_file)
        fits.writeto(small_file, small_image, useheader, overwrite=True)

        # plot_images(fits_path, data, small_image)

        print("    Reduction Complete!\n")
    print("```````````````````````````\n")

    return small_file

def plot_images(fits_path, data, small_image):
    import matplotlib.pyplot as plt
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


def load_magnetogram_params(datdir):
    """Reads the magnetic_target.params file and returns the parameters."""
    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'r') as fp:
        hdr = fp.readline().rstrip()
        cr = fp.readline().rstrip()
        fname = fp.readline().rstrip()
        adapt = int(fp.readline().rstrip())
        doplot = int(fp.readline().rstrip())
        reduce = int(fp.readline().rstrip())
    return (hdr, cr, fname, adapt, doplot, reduce)


def read_fits_data(fname):
    """Reads FITS data and fixes/ignores any non-standard FITS keywords."""
    hdulist = fits.open(fname, ignore_missing_simple=True)
    hdulist.verify('silentfix+warn')
    return hdulist

def load_fits_magnetogram(datdir = "/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/", batch="fluxon", bo=2, bn=2, ret_all=False):

    def read_fits_data(fname):
        """Reads FITS data and fixes/ignores any non-standard FITS keywords."""
        hdulist = fits.open(fname)
        hdulist.verify('silentfix+warn')
        return hdulist

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

def find_file_with_string(directory, search_string):
    for file_name in os.listdir(directory):
        if search_string in file_name:
            return os.path.join(directory, file_name)
    return None

def reduce_mag_file(mag_file, reduction=3, force=False):

    small_file = PosixPath(str(mag_file).replace("_r1_", f"_r{reduction}_"))
    # reduce the FITS image
    print(f"\tReducing image size by a factor of {reduction}...", end="")
    if not os.path.exists(small_file) or force:
        small_file = reduce_fits_image(mag_file, small_file, target_resolution=None, reduction_amount=reduction)
        # print("Success!\n")
    else:
        print("Skipped! Reduced file already exists:")
        print("\t\t", small_file)
        print("\n````````````````````````````````\n")

    return small_file

def write_params_file(datdir, cr, file_path, reduction):
    # write the parameter file
    params_path = os.path.join(datdir,"magnetic_target.params")
    with open(params_path, 'w') as fp:
        fp.write("## CR_int, Filename_str, Adapt_bool, Doplot_bool, reduction ##\n")
        fp.write(str(cr)+"\n")
        fp.write(str(file_path)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(0)+"\n")
        fp.write(str(reduction))



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
#     # mag_file = find_file_with_string(datdir, "br_r1")
#     # # reduce the FITS image
#     # small_file = reduce_mag_file(mag_file, reduction=reduction, force=force)
#     # # write the parameter file
#     # write_params_file(datdir, cr, small_file, reduction)
#     # # run the hilbert footpoint finder
#     # run_hilbert_footpoint_finder(datdir, cr, reduction=reduction, force=force)
#     # # return the hilbert footpoint file
#     # return find_file_with_string(datdir, "floc_cr"+str(cr))

#     # import os
#     # import numpy as np

#     # print("\n\n**Processing Magnetograms to get Footpoints...")

#     # no_floc = 0
#     # out_dir = os.path.join(datdir, batch_name, "cr" + str(cr), "floc")



# if __name__ == "__main__":
#     get_magnetogram_file(force_download=False, datdir="/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/")








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