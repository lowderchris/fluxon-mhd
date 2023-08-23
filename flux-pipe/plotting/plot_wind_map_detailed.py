"""_summary_

Raises:
    to: _description_

Returns:
    _type_: _description_
"""

import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
import matplotlib.pyplot as plt
import argparse
import os.path as path
from scipy.interpolate import griddata

from py_plot_helper import get_ax
from py_pipe_helper import \
    (load_fits_magnetogram, load_magnetogram_params, get_fixed_coords)
from plot_fieldmap import magnet_plot


# Remove outliers from the dataset recursively
def remove_outliers(data, ph, th, thresh_low=4, thresh_high=2, n_times= 1):
    """_summary_

    Args:
        data (_type_): _description_
        ph (_type_): _description_
        th (_type_): _description_
        thresh_low (int, optional): _description_. Defaults to 4.
        thresh_high (int, optional): _description_. Defaults to 2.
        n_times (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    # Get the Good Points
    mean = np.mean(data[data>0])
    std = np.std(data[data>0])
    data_good, inds_good = (list(t) for t in zip(*[(x,i) for i,x in enumerate(data)
                        if (mean - thresh_low * std < x < mean + thresh_high * std) and x > 0]))
    ph_good, th_good = ph[inds_good], th[inds_good]
    data_good = np.asarray(data_good)

    # Get the Bad Points
    inds_bad = [i for i in range(len(data)) if i not in inds_good]
    ph_bad, th_bad = ph[inds_bad], th[inds_bad]
    data_bad = data[inds_bad]

    # Return or Recurse
    if n_times <= 1:
        return data_good, ph_good, th_good, data_bad, ph_bad, th_bad
    else:
        data_good_2, ph_good_2, th_good_2, data_bad_2, ph_bad_2, th_bad_2 = \
                remove_outliers(data_good, ph_good, th_good, thresh_low, thresh_high, n_times-1)
        bad_data_all = np.concatenate((data_bad, data_bad_2))
        bad_ph_all   = np.concatenate((  ph_bad,   ph_bad_2))
        bad_th_all   = np.concatenate((  th_bad,   th_bad_2))
        return data_good_2, ph_good_2, th_good_2, bad_data_all, bad_ph_all, bad_th_all

def scale_data(vel0_clean, vel1_clean, outlier_V0, outlier_V1, scale=15**2, power=1):
    """Scale the data between 0 and 1, then raise to a power, then scale by a factor

    Args:
        vel0_clean (_type_): _description_
        vel1_clean (_type_): _description_
        outlier_V0 (_type_): _description_
        outlier_V1 (_type_): _description_
        scale (_type_, optional): _description_. Defaults to 15**2.
        power (int, optional): _description_. Defaults to 1.

    Raises:
        to: _description_

    Returns:
        _type_: _description_
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



def hist_plot(vel1_clean, ax=None, vmin=400, vmax=800, n_bins=20, do_print_top=True):
    """_summary_

    Args:
        vel1_clean (_type_): _description_
        ax (_type_, optional): _description_. Defaults to None.
        vmin (int, optional): _description_. Defaults to 400.
        vmax (int, optional): _description_. Defaults to 800.
        n_bins (int, optional): _description_. Defaults to 20.
        do_print_top (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    if do_print_top:
        print("\n\t\tMaking Histogram Plot...", end='')

    ## The Histogram Plot
    fig, hist_ax = get_ax(ax)
    mean1 = np.mean(vel1_clean)
    median1 = np.median(vel1_clean)
    std1  =  np.std(vel1_clean)
    hist_ax.hist(vel1_clean[vel1_clean>=0], bins=n_bins, color='sandybrown', density=True)
    hist_ax.axvline(mean1, color='k', linestyle='dashed', linewidth=1,
                    label=f"Mean: {mean1:.0f} km/s")
    hist_ax.axvline(median1, color='lightgrey', linestyle='-.', linewidth=1,
                    label=f"Median: {median1:.0f} km/s")
    hist_ax.axvline(mean1+std1, color='k', linestyle=':', linewidth=1,
                    label=f"Std: {std1:.0f} km/s")
    hist_ax.axvline(mean1-std1, color='k', linestyle=':', linewidth=1)
    hist_ax.legend()
    hist_ax.set_xlabel("Velocity (km/s)")
    hist_ax.set_ylabel("Fluxon Relative Frequency")
    fig.suptitle(f'CR{CR}, {len(vel1_clean)} Open Fields')

    hist_ax.set_xlim((vmin, vmax))
    hist_ax.set_ylim((0, 0.02))
    if do_print_top:
        print("Success!")
    return mean1, std1

def magnet_plot_orig(batch, ax=None, doplot=False, vmin=-500, vmax=500):
    """_summary_

    Args:
        batch (_type_): _description_
        ax (_type_, optional): _description_. Defaults to None.
        doplot (bool, optional): _description_. Defaults to False.
        vmin (int, optional): _description_. Defaults to -500.
        vmax (int, optional): _description_. Defaults to 500.

    Returns:
        _type_: _description_
    """
    magnet = load_fits_magnetogram(batch=batch, bo=3, bn=2)
    # magnet, header = load_fits_magnetogram(batch=batch, ret_all=True)
    # find the max and min of the magnetogram plot for use in setting the colormap,

    # Plot the magnetogram
    # Scatter the fluxon data
    ph000 = (ph0_clean+np.pi)%(2*np.pi)
    ph111 = (ph1_clean+np.pi)%(2*np.pi)
    ph000b = (ph0b+np.pi)%(2*np.pi)
    ph111b = (ph1b+np.pi)%(2*np.pi)
    sc00, sc01, magimg = None, None, None
    if doplot:
        _, ax = get_ax(ax)    ## The Magnetogram Plot
        magimg = ax.imshow(magnet, cmap='gray', interpolation=None,
                           origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto',
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

def hex_plot(ph1_clean, th1_clean, vel1_clean, ax=None, nx=20, vmin=400, vmax=800, do_print_top=True, do_hex=True):
    """_summary_

    Args:
        ph1_clean (_type_): _description_
        th1_clean (_type_): _description_
        vel1_clean (_type_): _description_
        ax (_type_, optional): _description_. Defaults to None.
        nx (int, optional): _description_. Defaults to 20.
        vmin (int, optional): _description_. Defaults to 400.
        vmax (int, optional): _description_. Defaults to 800.
        do_print_top (bool, optional): _description_. Defaults to True.
        do_hex (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    if do_print_top:
        if do_hex:
            print("\n\t\tMaking Hexbin Plot...", end="")
        else:
            print("\n\t\tMaking Contour Plot...", end="")
    _, hex_ax = get_ax(ax)
    ## Plot the Interp. data
    # Define the x-boundaries of the domain
    x_min = 0
    x_max = 2*np.pi

    # Wrap the data around the edges of the domain
    ph1 = (ph1_clean)%(2*np.pi)
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
    if do_hex:
        hex1 = hex_ax.hexbin(grid_x.flatten(), grid_y.flatten(), C=z1_use.flatten(),
                gridsize=gridsize, cmap='autumn',
                # vmin=np.nanmin(z1_use), vmax=np.nanmax(z1_use))
                vmin=vmin, vmax=vmax)
        print("Done!")
        return hex1
    else:
        # hex1 = hist_ax.imshow(grid_z1, extent=(0, 2*np.pi, -1, 1), zorder=0, alpha=1, cmap="autumn", vmin=vmin, vmax=vmax)
        # contour1 = hex_ax.contourf(grid_z1, zorder=-10, alpha=1, cmap="autumn")
        contour1 = hex_ax.contourf(grid_x, grid_y, grid_z1, zorder=0, alpha=1, cmap="autumn", vmin=vmin, vmax=vmax)
        print("Done!")
        return contour1


    if do_print_top:
        print("Success!")



if __name__ == "__main__":
## CODE STARTS HERE

    # create the argument parser
    parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
    parser.add_argument('--cr', type=int, default=None, help='Carrington Rotation')
    parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxons/fluxon-data', help='data directory')
    parser.add_argument('--show', type=int, default=1)
    parser.add_argument('--interp', type=str, default="linear")
    # parser.add_argument('--nact', type=int, default=0)
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
    vel0_clean, ph0_clean, th0_clean, v0b, ph0b, th0b = remove_outliers(vel0, ph0, th0, 3, 3, 1)
    vel1_clean, ph1_clean, th1_clean, v1b, ph1b, th1b = remove_outliers(vel1, ph1, th1, 3, 2, 3)
    v0, v1, v0bs, v1bs = scale_data(vel0_clean, vel1_clean, v0b, v1b, scale=15**2, power=1)

    ## PLOTTING
    fig, ax = plt.subplots(5)

    mag_ax = ax[0]
    scatter_ax = ax[1]
    contour_ax = ax[2]
    hex_ax = ax[3]
    hist_ax = ax[4]

    all_vmin, all_vmax = 475, 700
    drk=0.25
    n_open, n_closed, n_flux, fnum, n_outliers = magnet_plot(CR, dat_dir, batch, ax=mag_ax, vmin=-500, vmax=500, reduce_amt=reduce, nwant=args.nwant, do_print_top=False)
    hex_n = n_open//10
    hex1 = hex_plot(ph1_clean, th1_clean, vel1_clean, ax=hex_ax, nx=hex_n, vmin=all_vmin, vmax=all_vmax)

    scatter_ax.set_facecolor('grey')
    contour_ax.set_facecolor('grey')

    scat1 = scatter_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=6**2, alpha=0.75, marker='s', vmin=all_vmin, vmax=all_vmax, cmap='autumn', edgecolors='none')
    cont1 = contour_ax.scatter(ph1_clean, th1_clean, c=vel1_clean, s=4**2, alpha=1.0, marker='o', vmin=all_vmin, vmax=all_vmax, cmap='autumn', edgecolors='none')

    scatter_ax.scatter(ph1b, th1b,           color='g', s=3**2, alpha=1, marker='+')
    contour_ax.scatter(ph1b, th1b,           color='g', s=3**2, alpha=1, marker='+')
    scatter_ax.scatter(ph1_clean, th1_clean, c='k', s=2**2, alpha=1., marker='o')


    n_bins = np.linspace(all_vmin, all_vmax, 32)
    mean1, std1 = hist_plot(vel1_clean, ax=hist_ax, vmin=all_vmin, vmax=all_vmax, n_bins=n_bins)


    ## SAVING
    # Set the output file names
    filename = f"png_cr{CR}_f{args.nwant}_ou{n_open}_radial_wind.png"
    main_file =  f'{dat_dir}/batches/{batch}/cr{CR}/wind/{filename}'
    outer_file = f"{dat_dir}/batches/{batch}/imgs/windmap/{filename}"

    import os
    if not path.exists(os.path.dirname(main_file)):
        os.makedirs(os.path.dirname(main_file))

    if not os.path.exists(os.path.dirname(outer_file)):
        os.makedirs(os.path.dirname(outer_file))

    for this_ax in [mag_ax, hex_ax, scatter_ax, contour_ax]:
        this_ax.set_ylabel('sin(latitude)')
        this_ax.set_ylim((-1.0,1.0))
        this_ax.set_aspect('equal')
        this_ax.axhline(-1, c='lightgrey', zorder=-10)
        this_ax.axhline( 1, c='lightgrey', zorder=-10)
        this_ax.axvline(0, c='lightgrey', zorder=-10)
        this_ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
        this_ax.set_xlim((0, 2*np.pi))

    for jj in [0,1,2]:
        ax[jj].set_xticklabels([])

    # Add an axes for the colorbar
    cbar_ax = fig.add_axes([0.88, 0.2, 0.03, 0.75])

    # Create a colorbar with a custom colormap including the green overlay
    cmap = mpl.colormaps['autumn']

    plotobjs = [scat1, cont1, hex1]
    for obj in plotobjs:
        cbar = plt.colorbar(obj, cax=cbar_ax, extend="both", cmap=cmap, extendfrac=0.1,
                            aspect=15)
        cbar.cmap.set_over('lime')
        cbar.cmap.set_under('darkviolet')
        # cbar.cmap.set_over('lightgreen')
        # cbar.cmap.set_under('lightblue')

    cbar.set_label("Interp. Wind Speed [km/s]", labelpad=-50)

    fig.set_size_inches((6,8))

    plt.subplots_adjust(hspace=0.0,)

    # Save the Figures
    print("\n\t\tSaving figures to disk...", end="")
    # print(main_file)
    main_pdf = main_file.replace(".png", ".pdf")
    outer_pdf = outer_file.replace("png", "pdf")

    plt.savefig(outer_file, dpi=200)
    # print("\t\t\tSaving ", shorten_path(outer_pdf))
    plt.savefig(outer_pdf, dpi=200)

    plt.close(fig)
    print("Success!")

    print("\n\t    Done with wind plotting!\n")
    print("\t\t\t```````````````````````````````\n\n\n")
