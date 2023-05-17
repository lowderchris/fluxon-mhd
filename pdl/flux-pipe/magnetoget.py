# Import libraries
import sunpy
import sunpy.io
import sunpy.coordinates
from sunpy.net import Fido, attrs as a
import drms
import os
import glob



def get_magnetogram_files(cr=None, date=None, data_dir=None, email=None, do_download=True):
    """
    Function to grab and update MDI and HMI data.

    Args:
        cr (int): Carrington rotation number.
        date (str): Date in YYYY-MM-DD format.
        data_dir (str): Optional directory where data will be stored. If not specified,
            default directories will be used.

    Returns:
        None
    """


    if cr is not None:
        which = f"CR{cr}"
    elif date is not None:
        which = date
    else:
        which = 2193
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"\n\tDownloading Magnetograms for {which}...")

    # print("\nConfiguring Download...")
    try:
        jsoc_email = email or os.environ["JSOC_EMAIL"]
    except KeyError:
        jsoc_email = "chris.gilly@colorado.edu"

    # Specify any directories
    if data_dir is not None:
        if cr is not None:
            data_dir += f"/fluxon/cr{cr}/mag"
        hmidat = os.path.expanduser(os.path.join(data_dir, 'hmi.Synoptic_Mr.polfil'))
        mdidat = os.path.expanduser(os.path.join(data_dir, 'mdi.Synoptic_Mr.polfil'))
    else:
        hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil')
        mdidat = os.path.expanduser('~/data/mdi.Synoptic_Mr.polfil')

    # HMI data

    # Sort out the last downloaded rotation
    crfiles = glob.glob(hmidat+'*.fits')
    crfiles.sort()
    crlist = [int(i[-19:-15]) for i in crfiles]

    # Specify requested rotations
    if cr is None and date is None:
        if len(crlist) == 0:
            cr0 = 2096
        else:
            cr0 = max(crlist) + 1
        cr1 = int(sunpy.coordinates.sun.carrington_rotation_number(t='now'))

        if (cr0 - 1) == cr1:
            print('what?')
    else:
        if cr is not None:
            cr0 = cr
        else:
            cr0 = int(sunpy.coordinates.sun.carrington_rotation_number(date))

        cr1 = cr0

    # Start the client
    load_failed = False
    if not do_download:
        try:
            import pathlib
            use_path = which
            # print(which)
            if not "CR" in use_path:
                use_path = "CR" + use_path
            if not "mag" in data_dir:
                data_dir = data_dir + f"/fluxon/{use_path}/mag/"
            if not "hmi.Synoptic_Mr.polfil" in data_dir:
                data_dir = data_dir + "/hmi.Synoptic_Mr.polfil"
            hmi_object = pathlib.Path(data_dir)
            # print(hmi_object)
            file_list = list(hmi_object.iterdir())
            for file in file_list:
                ff = str(file)
                if 'fits' in ff and 'small' not in ff:
                    hmi_path_out = ff
        except FileNotFoundError as e:
            # print(e)
            # print("This file hasn't been downloaded yet!")
            load_failed = True

        if do_download or load_failed:
            print("\n\tDownloading HMI from JSOC...")
            c = drms.Client()

            # Generate a search
            crots = a.jsoc.PrimeKey('CAR_ROT', str(cr0) + '-' + str(cr1))
            res = Fido.search(a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), crots, 
                            a.jsoc.Notify(jsoc_email))

            # Once the query is made and trimmed down...
            hmi_path = hmidat+'/{file}.fits'
            hmi_path_out = Fido.fetch(res, path=hmi_path)[0]
    
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

    print("\tDownload Complete!\n")

    return (hmi_path_out, None)


from astropy.nddata import block_reduce
from astropy.io import fits
import numpy as np


def reduce_fits_image(fits_path, target_resolution=None, reduction_amount=None, func=np.nanmean):
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
    print(f"\tReducing {fits_path}...")
    # Open the FITS file and read the data
    with fits.open(fits_path, ignore_missing_simple=True) as hdul:
        hdul.verify('silentfix')
        data = hdul[0].data
        if data is None:
            data = hdul[1].data
            
        current_resolution = max(data.shape)
        print("\tOriginal Size: ", data.shape)

        # Calculate the reduction amount if target resolution is specified
        if target_resolution is not None:
            reduction_amount = int(np.ceil(current_resolution / target_resolution))

        # Raise an error if neither target_resolution nor reduction_amount is specified
        elif reduction_amount is None:
            raise ValueError("Either target_resolution or reduction_amount must be specified.")

        # Reduce the image and save it to a new file with "_small" appended to the filename
        small_image = block_reduce(data, reduction_amount, func)
        output_path = fits_path.replace('.fits', '_small.fits')
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

        fits.writeto(output_path, small_image, useheader, overwrite=True)

        print("\tFinal Size:    ", small_image.shape)
        # print("CDELTA 1 and 2: ", useheader['CDELT1'], useheader['CDELT2'])
        # print(small_image.shape[0] * useheader['CDELT2'] * np.pi/2)
        # print(useheader['CDELT1'] * small_image.shape[1])

        plot_images(fits_path, data, small_image)

        print("\n\tReduction Complete!\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    return output_path

def plot_images(fits_path, data, small_image):
    import matplotlib.pyplot as plt
    # Save the high resolution image as a grayscale PNG
    plt.imshow(data, cmap='gray')
    plt.axis('off')
    high_res_output_path = fits_path.replace('.fits', '.png')
    fig = plt.gcf()
    shp = data.shape

    ratio = shp[1]/shp[0]
    sz0=6 #inches
    sz1=sz0*ratio #inches
    DPI = shp[1] / sz1 #pixels/inch
    fig.set_size_inches((sz1, sz0))
    plt.savefig(high_res_output_path, bbox_inches='tight', dpi=4*DPI)
    plt.close()

    # Save the low resolution image as a grayscale PNG
    plt.imshow(small_image, cmap='gray')
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
    params_path = datdir + "magnetic_target.params"
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
    hdulist = fits.open(fname)
    hdulist.verify('silentfix+warn')
    return hdulist

def load_fits_magnetogram(datdir = "/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/"):

    def read_fits_data(fname):
        """Reads FITS data and fixes/ignores any non-standard FITS keywords."""
        hdulist = fits.open(fname)
        hdulist.verify('silentfix+warn')
        return hdulist

    fname = load_magnetogram_params(datdir)[2]
    fits_path = datdir + fname
    hdulist = read_fits_data(fits_path)
    brdat = hdulist[0].data
    header= hdulist[0].header
    brdat = brdat - np.mean(brdat)

    return brdat


if __name__ == "__main__":
    get_magnetogram_files(do_download=False, data_dir="/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/")
