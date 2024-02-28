from astropy.io import fits
from scipy.ndimage import zoom
import numpy as np
import matplotlib as mpl
mpl.use('qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from skimage.measure import block_reduce
import sunpy.sun.constants as const

import cv2
## Change the labels to show that the plots are inunits of std

root = "/Users/cgilbert/vscode/fluxon-data/comp_mags_png_5/"
# The good sigma is 0.3061


# 720 x 288
# image1_path = "/Users/cgilbert/vscode/fluxon-data/magnetograms/CR2193_r5_hmi.fits"

# Some Resolution
image1_path = "/Users/cgilbert/vscode/fluxon-data/magnetograms/CR2193_r2_hmi.fits"

# 180 x 360
image2_path = "/Users/cgilbert/vscode/fluxon-data/magnetograms/ADAPT/CR2193_rf2_adapt40311_03k012_201707200800_i00015600n1.fits"


def plot_scatter(image1, image2):
    plt.figure(figsize=(6,6))
    plt.scatter(image1.flatten(), image2.flatten(), alpha=0.1, edgecolors='w', linewidths=0.5)
    plt.xlabel('Image 1')
    plt.ylabel('Image 2')
    plt.title('Scatter plot of pixel intensities')
    # plt.show()

def downsample_image(image, block_size):
    return block_reduce(image, block_size=block_size, func=np.sum)





#############################################################################################

from scipy.ndimage import map_coordinates, geometric_transform


# def inverse_sine_warp(coords):
#     # The incoming coordinates are in the form (y, x). 
#     # We apply the arcsin transformation only to the y-coordinate.
#     out_coords = np.asarray([np.arcsin(coords[1]), coords[0]])
#     # print(out_coords, " : ", type(out_coords))
#     # return out_coords
#     return np.asarray([float(np.arcsin(coords[1])), float(coords[0])])


# def warp_image2(image):
#     # print(type(inverse_sine_warp))
#     return geometric_transform(image, np.sin)

# Call this function to warp your image.
# warped_image2 = warp_image(image2)

from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp

def reproject_fits(input_fits, hduID=0):
    # Load the FITS hdulist using astropy.fits
    hdu = fits.open(input_fits)[hduID]

    # Get the WCS info from the FITS file
    wcs_in = WCS(hdu.header)

    # Define new WCS with SIN projection for latitude
    wcs_out = wcs_in.deepcopy()
    wcs_out.wcs.ctype[1] = 'DEC--SIN'

    # Shape for the output array
    shape_out = hdu.data.shape

    # Reproject the image data to the new WCS
    array_out, footprint = reproject_interp(hdu, wcs_out, shape_out=shape_out)



    # # Reproject the image data to the new WCS
    # array_out, footprint = reproject_interp(hdu, wcs_out)

    # Return the reprojected data
    return array_out, footprint


import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

def reproject_fits(input_fits):
    # Load the FITS hdulist using astropy.fits
    hdu = fits.open(input_fits)[0]

    # Get the WCS info from the FITS file
    wcs_in = WCS(hdu.header, naxis=2)

    # Create a grid of coordinates from the WCS
    y, x = np.indices(hdu.data.shape)
    ra, dec = wcs_in.all_pix2world(x, y, 0)

    # Convert to astropy SkyCoord object
    coords = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')

    # Create new WCS with SIN projection for latitude
    wcs_out = wcs_in.deepcopy()
    wcs_out.wcs.ctype[1] = 'DEC--SIN'

    # Convert coordinates to new system
    new_ra, new_dec = wcs_out.all_world2pix(coords.ra.deg, coords.dec.deg, 0)

    # Interpolate data to new coordinates
    from scipy.interpolate import griddata
    array_out = griddata((new_ra.flatten(), new_dec.flatten()), hdu.data.flatten(), (x.flatten(), y.flatten()), method='cubic')

    return array_out.reshape(hdu.data.shape)





# # Use the function
# reprojected_data = reproject_fits(input_fits)



def warp_image(image):
    # Create a grid of the same shape as the image
    x, y = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    # orig_coords = np.array([x, y])

    xx = (x - 90) / 90

    # Apply the inverse transformation to the y-coordinates
    x_transformed = np.arcsin(xx)
    x_transformed= ((x_transformed/(3.1415926)*180)+90)

    # Stack the coordinates together
    coords = np.array([x_transformed, y])

    # Use map_coordinates to apply the transformation
    warped_image = map_coordinates(image, coords, order=1)

    return warped_image


def blur_image(img, sigma):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(img, sigma=sigma)


def do_all(kernel_width = 5, warp=True, title=None):
    # Call this function to warp your image.
    # warped_image2 = warp_image(image2)

    # Load the images
    img1_hdu = fits.getdata(image1_path, header=True)
    img2_hdu = fits.getdata(image2_path, header=True)
    # image2 = warp_image(image2)
    # image2 = reproject_fits(image2_path)

    # std2 = np.std(image2)
    # mean2 = np.mean(image2)
    # vmin = -10 #mean2 - 0.1*std2
    # vmax = 10 #mean2 + 0.1*std2

    # plt.imshow(image2, cmap='hot', interpolation='nearest', origin='lower',
    #            vmin = vmin, vmax = vmax)
    # plt.show()


    # image2, header2 = fits.getdata(image2_path,  header=True)



    compareplot(img1_hdu, img2_hdu, kernel_width, title=title)


    # # Perform difference analysis
    # # difference = normalized_image1 - normalized_image2


    # # Check if the total sums are equal
    # print("Total sum comparison:", np.sum(downsampled_image1) == np.sum(image2))

    # print("Difference:", np.sum(np.abs(difference)))

    # # Additional analysis: Calculate mean and standard deviation
    # print("Mean comparison:", np.mean(normalized_image1) == np.mean(normalized_image2))

    # print("Standard Deviation comparison:", np.std(normalized_image1) == np.std(normalized_image2))

    # # # Visualize difference
    # # plt.imshow(difference, cmap='hot', interpolation='nearest')
    # # plt.colorbar()
    # # plt.show()

# import numpy as np
# import sunpy.constants as const

def calculate_pixel_areas(image_shape, long_range=(-180, 180), lat_range=(-90, 90)):
    # Calculate the differential in longitude and latitude (in radians)
    d_long = np.deg2rad(np.abs(long_range[1] - long_range[0]) / image_shape[1])
    d_lat = np.deg2rad(np.abs(lat_range[1] - lat_range[0]) / image_shape[0])
    
    # Calculate the latitude at the top and bottom of each pixel (in radians)
    latitudes_rad = np.linspace(np.deg2rad(lat_range[0]) + d_lat / 2, np.deg2rad(lat_range[1]) - d_lat / 2, image_shape[0])
    
    # Area of each pixel (differential area element on a sphere)
    pixel_areas = (const.radius**2 * d_long * np.abs(np.sin(latitudes_rad[:, None] + d_lat / 2) - np.sin(latitudes_rad[:, None] - d_lat / 2))).to('m^2')

    return pixel_areas

# # Example usage:
# image_shape = (1024, 1024)  # Replace with your actual image shape
# pixel_areas = calculate_pixel_areas(image_shape)

# print(pixel_areas)


def compareplot(img1_hdu, img2_hdu, width, warp=True, title=None):
    image1, header1 = img1_hdu
    image2, header2 = img2_hdu


    # Make image1 be in units of Maxwells
    surface_area_sun = const.surface_area.to(u.cm*u.cm)
    n_pix = image1.shape[0] * image1.shape[1] 
    area_per_pix = surface_area_sun / n_pix
    image1 = image1 * u.Maxwell / u.cm**2
    image1 *= area_per_pix



    #### Make image2 be in units of Maxwells
    n_pix2 = image2.shape[0] * image2.shape[1] 

    # find the area per pixel for an image that is sin(latitude) in y and (longitude) in x
    # This is the same as the area per pixel for the warped image
    # The result should be an array with each element equal to the area for that pixel

    pixel_areas = calculate_pixel_areas(image2.shape, long_range=(0, 360), lat_range=(0, 180)) #m**2
    pixel_trans_areas = np.zeros((180, 360))    # convert the 1d array to a 2d array with every column the same
    for i in np.arange(360):
        pixel_trans_areas[:, i] = pixel_areas[:, 0] # m**2

    # surfaceArea_pixels = np.sum(pixel_trans_areas) * u.m**2

    united_im2 = (image2 * u.Maxwell/u.cm**2) * (pixel_trans_areas * u.m**2).to(u.cm**2) # Maxwells

    # pixel_areas_transposed = np.tile(pixel_areas, (image2.shape)).T



    # Multiply the image by the area per pixel to get the total flux in each pixel  
    if warp:
        warped_united_im2 = warp_image(united_im2)

    # thing = np.sum(np.abs(warped_united_im2), axis=1)
    # unwarped_thing = np.sum(np.abs(united_im2), axis=1)
    # plt.plot(thing)
    # plt.plot(unwarped_thing)
    # plt.show()
    # plt.imshow(warped_united_im2, origin='lower')
    # plt.show()

    # kernel_width = 10
    # blurred_img1 = cv2.GaussianBlur(image1, (kernel_width, kernel_width), 0)
    blurred_img1 = blur_image(image1, sigma=width)

    # Get the ratio of sizes
    ratio = np.array(image1.shape) / np.array(image2.shape)

    # Downsample the larger image to the size of the smaller one
    downsampled_image1 = zoom(blurred_img1, 1 / ratio)

    asdf = 1
    # print(downsample_image.shape)
    # compareplot(downsampled_image1, image2)


    # Normalize both images
    normalized_image1 = downsampled_image1  / np.nanmean(np.abs(downsampled_image1))
    normalized_image2 = image2              / np.nanmean(np.abs(image2))
    from copy import copy

    # from scipy.ndimage import gaussian_filter
    # image1 = gaussian_filter(image1, sigma=2)
    img_1_orig = copy(image1)

    std1 = np.std(image1)
    mean1 = np.mean(image1)
    std2 = np.std(image2)
    mean2 = np.mean(image2)
    image1 -= mean1
    image1 /= 2*std1
    image2 -= mean2
    image2 /= std2
    difference = image1-image2
    signed_diff = np.sum(difference)
    unsigned_diff = np.sum(np.abs(difference))
    ### PLOTTING ###

    # Create 4 subplots
    gs = gridspec.GridSpec(4, 1)  # Adjust ratios as necessary
    fig = plt.figure(figsize=(7, 14))
    ax = [plt.subplot(gs[0,0])]
    for i in np.arange(1, 3):
        ax.append(plt.subplot(gs[i], sharex=ax[0], sharey=ax[0]))
    ax_hist = plt.subplot(gs[3])
    ax.append(ax_hist)


    # Plot the images
    vmax = np.mean((std1, std2))
    vmin = -np.mean((std1, std2))
    im0 = ax[0].imshow(image1, cmap='hot', interpolation='nearest', vmin=vmin, vmax=vmax)
    im1 = ax[1].imshow(image2, cmap='hot', interpolation='nearest', vmin=vmin, vmax=vmax)
    im2 = ax[2].imshow(difference, cmap='hot', interpolation='nearest', vmin=vmin, vmax=vmax)


    # Call this function after you have normalized your images
    # plot_scatter(image1, image2)





    # Flatten the images to create 1D arrays
    flat_normalized_image1 = image1.flatten()
    flat_normalized_image2 = image2.flatten()
    flat_normalized_difference = flat_normalized_image1 - flat_normalized_image2

    ## Plot histograms
    range = (-3, 3)
    bins = 50

    # Create Histograms
    hist_data1, bin_edges1 = np.histogram(flat_normalized_image1, bins=bins, range=range)  # Replace "data" with your histogram data
    hist_data2, bin_edges2 = np.histogram(flat_normalized_image2, bins=bins, range=range)  # Replace "data" with your histogram data
    hist_data_diff, bin_edges_diff = np.histogram(flat_normalized_difference, bins=bins, range=range)  # Replace "data" with your histogram data

    # Normalize the histograms
    normfunc = np.max
    hist_data1 = hist_data1 / normfunc(hist_data1)
    hist_data2 = hist_data2 / normfunc(hist_data2)
    hist_data_diff = hist_data_diff / normfunc(hist_data_diff)

    # hist_curves_diff = np.sqrt(np.abs(hist_data1**2 - hist_data2**2))
    hist_curves_diff = np.abs(hist_data1 - hist_data2)

    def smooth_array(array, window_size):
        window = np.ones(window_size) / window_size
        smoothed_array = np.convolve(array, window, mode='same')
        return smoothed_array

    hist_curves_diff = smooth_array(hist_curves_diff, 3)
    hist_curves_sum = np.nansum(hist_curves_diff)

    # Get bin centers
    bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2
    bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2
    bin_centers_diff = (bin_edges_diff[:-1] + bin_edges_diff[1:]) / 2

    # Plot histogram as a line
    ax_hist.plot(bin_centers1, hist_data1, alpha=0.75, label='HMI', lw=4, c='C0')
    ax_hist.plot(bin_centers2, hist_data2, alpha=0.75, label='ADAPT', lw=4,  c='C1')
    ax_hist.plot(bin_centers_diff, hist_data_diff, alpha=0.75, label='Difference', lw=4, c='C2')

    ax_hist.plot(bin_centers2, hist_curves_diff, alpha=0.75, label='Hist Diff', lw=4,  c='C3')


    # Format the histogram plot
    ax_hist.set_title('Normalized Histogram of Image Values: Sum = {:0.5}'.format(hist_curves_sum))
    ax_hist.set_xlabel('Pixel value (sigma = {:0.4}, {:0.4})'.format(2*std1, std2))
    ax_hist.set_ylabel('Frequency')
    ax_hist.set_yscale('log')
    ax_hist.legend(loc='upper right')
    ax_hist.grid(True)
    ax_hist.set_xlim(*range)
    ax_hist.set_ylim((10**-5, 2))
    # plt.show()

    # Format the image plots
    fig.suptitle("Magnetogram Comparison")
    ax[0].set_title("HMI: CR2193_r2_hmi.fits : sigma={:0.4}, normsum={:0.4}, sum={:0.4}".format(float(width), np.sum(image1), np.sum(img_1_orig)))
    ax[1].set_title("ADAPT: CR2193_rf2_a....fits : sum={:0.4}".format(np.sum(image2)))
    ax[2].set_title(f"Difference: {signed_diff:0.03} signed, {unsigned_diff:0.03} unsigned")


    # Add an axes for the colorbar
    cbar_ax = fig.add_axes([0.87, 0.2, 0.03, 0.75])
    fig.subplots_adjust(top=0.96,
                        bottom=0.045,
                        left=0.095,
                        right=0.85,
                        hspace=0.1,
                        wspace=0.2)

    # Create a colorbar with a custom colormap including the green overlay
    import matplotlib as mpl
    cmap = mpl.colormaps['autumn']

    plotobjs = [im0, im1, im2]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend="both", cmap=cmap, extendfrac=0.1,
                            aspect=15)
        cbar.cmap.set_over('lime')
        cbar.cmap.set_under('darkviolet')
        # cbar.cmap.set_over('lightgreen')
        # cbar.cmap.set_under('lightblue')

    cbar.set_label("Magnetic Field [Normed Gauss]") #, labelpad=-50)



    # Show the Plot
    # fig.tight_layout()
    import time
    import os
    now = time.time()
    if not os.path.exists(root):
        os.makedirs(root)

    a=1

    if type(width) is int:
        savestring = f"{root}Compare_Magnetograms_{width}.png"
    else:
        wid_string = str("{:06.03f}".format(width)).replace('.', '_')
        savestring = f"{root}Compare_Magnetograms_{wid_string}.png"
    # plt.show()
    if title is not None:
        savestring = savestring.replace(".png", f"_{title}.png")

    show = False
    if show:
        fig.subplots_adjust(
                    top = 0.930,
                    hspace=0.310,
                    wspace=0.2)
        plt.show()
    else:
        plt.savefig(savestring, dpi=300)
        # print(f"Image {width} saved!")
    plt.close()


#############################################################################################
import cv2
import os

def images_to_video(input_dir, output_file, fps=30):
    # Get the list of image file names in the directory
    image_files = sorted(os.listdir(input_dir))

    # Get the first image to retrieve dimensions for the video

    png_images = [x for x in image_files if x.endswith(".png")]
    # print(png_images)
    first_image = cv2.imread(os.path.join(input_dir, png_images[0]))
    height, width, channels = first_image.shape

    # Initialize the video writer
    # fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Specify the video codec (e.g., 'mp4v', 'x264')
    fourcc = cv2.VideoWriter_fourcc(*'avc1')  # Specify the video codec (e.g., 'mp4v', 'x264')
    out_path = os.path.join(input_dir, output_file)
    video_writer = cv2.VideoWriter(out_path, fourcc, fps, (width, height))

    from tqdm import tqdm
    # Iterate over the image files and write them to the video
    for image_file in tqdm(png_images):
        image_path = os.path.join(input_dir, image_file)
        image = cv2.imread(image_path)
        video_writer.write(image)

    # Release the video writer
    video_writer.release()
    print(f"Video saved to {out_path}")


#############################################################################################
import numpy as np
from tqdm import tqdm

input_directory = root

if os.path.exists(input_directory):
    import shutil
    shutil.rmtree(input_directory, ignore_errors=True)
    # os.rmdir(input_directory, recursive=True)
os.makedirs(input_directory)

#############################################################################################
#### Run the Program ########################################################################
#############################################################################################

make_images = True
if make_images:
    for warp in [True]:
        sigma_list = [] #[0.4737]
        sigma_list.extend(np.linspace(0.0, 1.0, 50, endpoint=True))
        # print(sigma_list)
        for sigma in tqdm(sigma_list):
            do_all(sigma, warp=warp, title="warped" if warp else None)

do_video = True
if do_video:
    fps = 10
    output_video = f'output_{fps}.mp4'

    # if os.path.exists(os.path.join(input_directory, output_video)):
    #     os.remove(os.path.join(input_directory, output_video)

    images_to_video(input_directory, output_video, fps)