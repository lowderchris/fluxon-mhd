"""_summary_

Returns
-------
_type_
    _description_

Raises
------
to
    _description_
"""

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
import os.path as path
from scipy.interpolate import griddata

from py_plot_helper import get_ax
from py_pipe_helper import load_fits_magnetogram, load_magnetogram_params, get_fixed_coords
from plot_fieldmap import magnet_plot



def scale_data(vel0_clean, vel1_clean, outlier_V0, outlier_V1, scale=15**2, power=1):
    """ Scale the data between 0 and 1, then raise to a power, then scale by a factor

    Parameters
    ----------
    vel0_clean : _type_
        _description_
    vel1_clean : _type_
        _description_
    outlier_V0 : _type_
        _description_
    outlier_V1 : _type_
        _description_
    scale : _type_, optional
        _description_, by default 15**2
    power : int, optional
        _description_, by default 1

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    to
        _description_
    """

    vel0_max = np.nanmax(vel0_clean)
    vel1_max = np.nanmax(vel1_clean)
    vel0_min = np.nanmin(vel0_clean)
    vel1_min = np.nanmin(vel1_clean)

    v0 = scale * ((np.abs(vel0_clean) - vel0_min) / (vel0_max-vel0_min))**power
    v1 = scale * ((np.abs(vel1_clean) - vel1_min) / (vel1_max-vel1_min))**power

    outlier_V0_scaled = scale * ((np.abs(outlier_V0) - vel0_min) / (vel0_max-vel0_min))**power
    outlier_V1_scaled = scale * ((np.abs(outlier_V1) - vel1_min) / (vel1_max-vel1_min))**power

    return v0, v1, outlier_V0_scaled, outlier_V1_scaled

# Remove outliers from the dataset
def remove_outliers(data, ph0_temp, th0_temp, threshold=3):
    """_summary_

    Parameters
    ----------
    data : _type_
        _description_
    ph0_temp : _type_
        _description_
    th0_temp : _type_
        _description_
    threshold : int, optional
        _description_, by default 3

    Returns
    -------
    _type_
        _description_
    """

    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    filtered_data, good_points = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data)
                                    if mean - threshold * std < x < mean + threshold * std]))
    ph0c, th0c = ph0_temp[good_points], th0_temp[good_points]
    bad_points = [i for i in range(len(data)) if i not in good_points]
    ph0_b, th0_b = ph0_temp[bad_points], th0_temp[bad_points]
    outlier_data = data[bad_points]
    mean = np.mean(filtered_data)
    std = np.std(filtered_data)
    return np.asarray(filtered_data), mean, std, ph0c, th0c, outlier_data, ph0_b, th0_b


def hist_plot(vel1_clean, ax=None, vmin=400, vmax=800, n_bins=20, do_print_top=True):
    """_summary_

    Parameters
    ----------
    vel1_clean : _type_
        _description_
    ax : _type_, optional
        _description_, by default None
    vmin : int, optional
        _description_, by default 400
    vmax : int, optional
        _description_, by default 800
    n_bins : int, optional
        _description_, by default 20
    do_print_top : bool, optional
        _description_, by default True

    Returns
    -------
    _type_
        _description_
    """
    if do_print_top:
        print("\n\t\tMaking Histogram Plot...", end='')

    ## The Histogram Plot
    _, hist_ax = get_ax(ax)
    mean1 = np.mean(vel1_clean)
    median1 = np.median(vel1_clean)
    std1  =  np.std(vel1_clean)
    hist_ax.hist(vel1_clean[vel1_clean>=0], bins=n_bins, color='sandybrown')
    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1, label="Mean: {mean1:.0f} km/s")
    hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1, label="Median: {median1:.0f} km/s")
    hist_ax.axvline(mean1+std1, color='k', linestyle=':', linewidth=1, label="Std: {std1:.0f} km/s")
    hist_ax.axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
    hist_ax.legend()
    hist_ax.set_xlabel("Velocity (km/s)")
    hist_ax.set_ylabel("Number of Fluxons")
    hist_ax.set_title(f'CR{CR}, {len(vel1_clean)} Open Fields')

    hist_ax.set_xlim((vmin, vmax))
    if do_print_top:
        print("Success!")
    return mean1, std1
    # ax[3].set_title(f"Solar Wind Velocity Distribution: CR{CR}, Open Fluxons: {len(vel1)}", fontsize=16)

def magnet_plot_orig(batch, ax=None, doplot=False, vmin=-500, vmax=500):
    """_summary_

    Parameters
    ----------
    batch : _type_
        _description_
    ax : _type_, optional
        _description_, by default None
    doplot : bool, optional
        _description_, by default False
    vmin : int, optional
        _description_, by default -500
    vmax : int, optional
        _description_, by default 500

    Returns
    -------
    _type_
        _description_
    """
    magnet = load_fits_magnetogram(batch=batch, bo=3, bn=2)
    # magnet, header = load_fits_magnetogram(batch=batch, ret_all=True)

    # Plot the magnetogram
    # Scatter the fluxon data
    ph000 = (ph0_clean+np.pi)%(2*np.pi)
    ph111 = (ph1_clean+np.pi)%(2*np.pi)
    ph000b = (ph0b+np.pi)%(2*np.pi)
    ph111b = (ph1b+np.pi)%(2*np.pi)
    sc00, sc01, magimg = None, None, None
    if doplot:
        _, ax = get_ax(ax)    ## The Magnetogram Plot
        magimg = ax.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto',
                        #    vmin=mvmin, vmax=mvmax)
                           vmin=vmin, vmax=vmax)
        sc00 = ax.scatter(ph000, th0_clean, c= vel0_clean, s=v0, alpha=0.5, label=r'V(1.0R$_\odot$)'  , cmap="winter",  )
        sc01 = ax.scatter(ph111, th1_clean, c= vel1_clean, s=v1, alpha=0.5, label=r'V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
        # sc01 = ax[1].scatter(ph1, th1, c= vel1, s=v1, alpha=0.75, label='V(21.5R$_\odot$)' , cmap="autumn",  marker='s')
        # sc10 = ax[1].scatter(ph0, th0, c= fr0 , s=f0, alpha=0.75, label='Fr(1.0R$_\odot$)' , cmap="winter", )
        # sc11 = ax[1].scatter(ph1, th1, c= fr1 , s=f1, alpha=0.75, label='Fr(21.5R$_\odot$)', cmap="autumn", marker='s')

        #Scatter the Outliers
        ax.scatter(ph000b, th0b, c= v0bs, alpha=0.6, label=r'V(1.0R$_\odot$)', cmap="winter", marker='X',  s=50, edgecolors='k')
        ax.scatter(ph111b, th1b, c= v1bs, alpha=0.6, label=r'V(1.5R$_\odot$)', cmap="autumn", marker='X',  s=50, edgecolors='k')
        # ax[1].scatter(ph000b, th0b, c= outlier_V0_scaled, alpha=0.95, label='V(1.0R$_\odot$)', cmap="winter", marker='X', s=50, edgecolors='k')
        # ax[1].scatter(ph111b, th1b, c= v1bs, alpha=0.95, label='V(1.5R$_\odot$)', cmap="autumn", marker='X', s=50, edgecolors='k')
    return ph000, ph000b, ph111, ph111b, sc00, sc01, magimg

def hex_plot(ph1_clean, th1_clean, vel1_clean, ax=None, nx=20, vmin=400, vmax=800, do_print_top=True):
    """_summary_

    Parameters
    ----------
    ph1_clean : _type_
        _description_
    th1_clean : _type_
        _description_
    vel1_clean : _type_
        _description_
    ax : _type_, optional
        _description_, by default None
    nx : int, optional
        _description_, by default 20
    vmin : int, optional
        _description_, by default 400
    vmax : int, optional
        _description_, by default 800
    do_print_top : bool, optional
        _description_, by default True
    """
    if do_print_top: print("\n\t\tMaking Hexbin Plot...", end="")
    fig, hex_ax = get_ax(ax)
    ## Plot the Interp. data
    # Define the x-boundaries of the domain
    x_min = 0
    x_max = 2*np.pi

    # Wrap the data around the edges of the domain
    ph1 = (ph1_clean+np.pi)%(2*np.pi)
    th1 = th1_clean
    vel1 = vel1_clean

    ph1_wrapped = np.concatenate((ph1 - x_max, ph1, ph1 + x_max))
    th1_wrapped = np.concatenate((th1, th1, th1))
    vel1_wrapped = np.concatenate((vel1, vel1, vel1))

    # Create a grid for interpolation
    Ny, Nx = load_fits_magnetogram(batch=batch, bo=3, bn=2).shape
    grid_x, grid_y = np.linspace(x_min, x_max, Nx, endpoint=False), np.linspace(-1, 1, Ny)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)

    # Interpolate values on the grid
    points1 = [(ph, th) for ph, th in zip(ph1_wrapped, th1_wrapped)]
    grid_z1 = griddata(points1, vel1_wrapped, (grid_x, grid_y), method=interp, fill_value=0)

    # Plot the interpolated data
    sy, sx = grid_x.shape
    ratio = float(sy)/float(sx)
    gridsize = (int(nx*1.2), int(np.round(nx*ratio)))
    grid_z1_NANed = grid_z1.copy()
    grid_z1_NANed[grid_z1_NANed==0] = np.nan

    z1_use = grid_z1_NANed
    hex1 = hex_ax.hexbin(grid_x.flatten(), grid_y.flatten(), C=z1_use.flatten(),
            gridsize=gridsize, cmap='autumn',
            vmin=vmin, vmax=vmax)

    # Add an axes for the colorbar
    cbar_ax = fig.add_axes([0.88, 0.38, 0.01, 0.3])

    # Add a colorbar to the figure
    cbar = fig.colorbar(hex1, cax=cbar_ax, extend="max")
    cbar.cmap.set_over('lightgreen')
    cbar.set_label("Interp. Wind Speed [km/s]")

    if do_print_top: print("Success!")

## CODE STARTS HERE

if __name__ == "__main__":

    # create the argument parser
    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/fluxons/fluxon-data', help='data directory')
    parser.add_argument('--show', type=int, default=1)
    parser.add_argument('--interp', type=str, default="linear")
    parser.add_argument('--nact', type=int, default=0)
    parser.add_argument('--nwant', type=int, default=0)
    parser.add_argument('--file', type=str, default=None, help='select the file name')
    parser.add_argument('--batch', type=str, default="fluxon_paperfigs_5", help='select the batch name')
    args = parser.parse_args()
    batch = args.batch
    interp = args.interp
    dat_dir = args.dat_dir

    # Load the magnetogram parameters
    (hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(dat_dir)
    CR = args.cr or cr


    # Load the wind file
    dat_file = args.file or f'{dat_dir}/batches/{batch}/cr{CR}/wind/cr{CR}_f{args.nwant}_radial_wind.dat'
    # print("loading file: ", dat_file)
    arr = np.loadtxt(dat_file).T
    try:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1, fr1b = arr
    except ValueError:
        fid, phi0, theta0, phi1, theta1, vel0, vel1, fr0, fr1 = arr
    try:
        nfluxon = arr.shape[1]
    except IndexError as e:
        print("IndexError: ", e)
        print("Wind Calculation Failed, Rerun.")
        exit()


    # Convert coords to correct coords
    ph0, th0 = get_fixed_coords(phi0, theta0)
    ph1, th1 = get_fixed_coords(phi1, theta1)


    print("\n\tPlotting Windmap...", end="\n" if __name__=="__main__" else "")

    # Get the Data
    vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 2, 1)
    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)
    v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

    ## PLOTTING
    fig, ax = plt.subplots(3)

    mag_ax = ax[0]
    hex_ax = ax[1]
    hist_ax = ax[2]

    all_vmin, all_vmax = 450, 850
    drk=0.25

    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(CR, dat_dir, batch,
        ax=mag_ax, vmin=-500, vmax=500, reduce_amt=reduce, nwant=args.nwant, do_print_top=False)
    hex_plot(ph1_clean, th1_clean, vel1_clean, ax=hex_ax, nx=20, vmin=all_vmin, vmax=all_vmax)
    mean1, std1 = hist_plot(vel1_clean, ax=hist_ax, vmin=all_vmin, vmax=all_vmax, n_bins=16)

    sc01 = mag_ax.scatter(ph1_clean, th1_clean, c='k', s=2**2, alpha=0.6, marker='+')
    sc01 = mag_ax.scatter(ph1b, th1b, color=(drk, drk, drk), s=2**2, alpha=0.6, marker='+')


    ## SAVING
    # Set the output file names
    filename = f"cr{CR}_f{args.nwant}_ou{n_open}_radial_wind.png"
    main_file =  f'{dat_dir}/batches/{batch}/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    import os
    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))



    for this_ax in [mag_ax, hex_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.set_ylim((-1.0,1.0))
        this_ax.set_aspect('equal')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline( 1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        this_ax.set_xlim((0, 2*np.pi))

    ax[0].set_xlabel('Longitude (Radians)')
    ax[1].set_xticklabels([])


    fig.set_size_inches((6,6))
    # plt.tight_layout()

    plt.subplots_adjust(
    top=0.981,
    bottom=0.097,
    left=0.113,
    right=0.869,
    hspace=0.312,
    wspace=0.18
    )

    # plt.show()

    # Save the Figures
    print("\n\t\tSaving figures to disk...")
    # print(main_file)
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace(".png", ".pdf")

    plt.savefig(outer_pdf, dpi=200)

    plt.close(fig)
    print("Success!")

    print("\n\t    Done with wind plotting!\n")
    print("\t\t\t```````````````````````````````\n\n\n")
