import numpy as np
from scipy import interpolate
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import os
import re
np.seterr(divide='raise',invalid='raise')



def wind_radial_onefluid(nr=90000,rxmax=2150,mu=0.5,Tmax=1.0e6,psi=4.2,delta=0.0):
    xRsun = 6.96e10
    boltzk = 1.380622e-16
    xmHyd = 1.677333e-24
    Gconst = 6.6732e-8
    xMsun = 1.989e33

    # set up radial arrays ("rx" is in units of solar radii)

    rx = 10.**(np.linspace(np.log10(1.000001),np.log10(rxmax),num=nr))
    rcm = rx * xRsun
    drcm_mid = np.gradient(rcm)

    # use the Cranmer & Schiff (2016) parameters to set up a T(r) function

    xstar = (psi-2.)**(1./(psi-3.))
    zeta0 = (xstar**(psi-2.))/(psi-3.)
    T = Tmax * (zeta0*(1./rx - rx**(2.-psi)))**delta

    # now that we know T(r), set up the sound speed and right-hand side of the Parker equation of motion

    a2 = boltzk*T/mu/xmHyd
    da2dr = np.gradient(a2)/drcm_mid
    RHS = -da2dr + 2.*a2/rcm - Gconst*xMsun/rcm/rcm
    dRHSdr = np.gradient(RHS)/drcm_mid

    # integrate the RHS to obtain the Kopp & Holzer (1976) quantity that is minimized at the critical point

    RHSINT = np.cumsum(RHS*drcm_mid)
    ircrit = np.argmin(RHSINT)
    interp = interpolate.interp1d(RHS[ircrit-5:ircrit+5],rx[ircrit-5:ircrit+5],kind="linear")
    rxcrit = interp(0.)

    # at the critical point, interpolate to find quantities needed for L'Hopital's rule

    interp = interpolate.interp1d(rx,a2,kind="linear")
    a2crit = interp(rxcrit)
    interp = interpolate.interp1d(rx,da2dr,kind="linear")
    da2dr_crit = interp(rxcrit)
    interp = interpolate.interp1d(rx,dRHSdr,kind="linear")
    dRHSdr_crit = interp(rxcrit)

    # apply L'Hopital's rule to find the slope du/dr at the critical point

    ucrit = np.sqrt(a2crit)
    ducdr = da2dr_crit/2./ucrit
    dudr_crit = 0.5*(ducdr + np.sqrt(ducdr**2 + 2.*dRHSdr_crit))

    # find the array indices that surround the critical point

    indr = np.linspace(0,nr-1,num=nr,dtype=int)
    indcritLO = np.max(indr.compress((rx < rxcrit).flat))
    indcritHI = indcritLO + 1
    drLO = xRsun * (rxcrit - rx[indcritLO])
    drHI = xRsun * (rx[indcritHI] - rxcrit)

    # set up the wind speed array u(r), and get ready to integrate up and down from critical point

    u = np.zeros(nr)
    u[indcritLO] = ucrit - drLO*dudr_crit
    u[indcritHI] = ucrit + drHI*dudr_crit

    # integrate up!

    for i in range(indcritHI,nr,1):
        drcm = rcm[i]-rcm[i-1]
        DDD = u[i-1] - a2[i-1]/u[i-1]
        u[i] = u[i-1] + drcm*RHS[i-1]/DDD

    # integrate down!

    for i in range(indcritLO,-1,-1):
        drcm = rcm[i+1]-rcm[i]
        DDD = u[i+1] - a2[i+1]/u[i+1]
        u[i] = u[i+1] - drcm*RHS[i+1]/DDD

    return rxcrit,rx,T,u,ucrit


















# ########################################################################################
# ########################################################################################
# ########################################################################################
# ########################################################################################


def load_data(filename, print_data=False):
    data = np.load(filename)
    if print_data:
        a = [print(d, '\t', data[d]) for d in data.keys()]
        print()
    return data

def plot_data(xx, mydict):
    for key in mydict.keys():
        yy = np.asarray(np.squeeze(mydict[key]))
        try:
            plt.plot(xx, yy, label=key)
        except ValueError:
            # print("Error with key: ", key)
            # plt.axvline(yy)
            pass
    plt.legend()
    plt.yscale('symlog')
    plt.xscale('log')
    plt.show()

def load_tempest(file):
    file= file.replace("tempest.dat", "tempest_result")

    # directory = configs['data_dir']
    # file = "T3"
    inputs = f"{file}_inputs.npz"
    miranda =f"{file}_miranda.npz"
    prospero=f"{file}_prospero.npz"
    fullRHS =f"{file}_fullRHS.npz"

    in_dict = load_data(inputs)
    mir_dict = load_data(miranda)
    pro_dict = load_data(prospero)
    rhs_dict = load_data(fullRHS)

    zx = in_dict['zx']
    zTR = in_dict['zTR']
    zcrit = mir_dict['zcrit']
    u_miranda = np.squeeze(mir_dict['u']) / 100 / 1000
    u_prospero = np.squeeze(pro_dict['u']) / 100 / 1000
    rho = np.squeeze(rhs_dict['rho'])
    return zx, zTR, zcrit, u_miranda, u_prospero, rho

def plot_tempest(file, ax=None):
    print(f"Plotting {file}")
    this_dir = file.split('/data/')[0]+'/imgs/wind/tempest/'
    if not os.path.exists(this_dir):
        os.makedirs(this_dir)

    # Define the regex pattern to find 'cr' followed by numbers
    pattern = r'cr(\d+)'

    # Search for the pattern in the given file path
    cr = re.search(pattern, file).group(1)

    zx, zTR, zcrit, u_miranda, u_prospero, rho = load_tempest(file)

    if ax is None:
        fig, ax = plt.subplots()

    # ax.plot(zx-1, u_miranda.T, 'g')
    # ax.plot(zx-1, u_miranda[0].T, 'g', label='$u_{sw}~$Miranda')
    ax.plot(zx, u_prospero.T, 'm:')
    ax.plot(zx, np.ones_like(u_prospero[0].T), 'm:', label='$u_{sw}~$Prospero')
    # ax.plot(zx, rho.T*10**5, label='rho x $10^5$')
    # print(zTR)
    Rsun = 6.955e10
    tr_label = f'TR Mean={np.mean(zTR):0.4}'
    crit_label = f'CR Mean={(np.mean(zcrit)/Rsun):0.4f}'
    for TR in zTR:
        ax.axvline(TR, linestyle='--', label=tr_label, zorder=-10)
        tr_label = None
    for crit in zcrit:
        ax.axvline(crit/Rsun-1, color='r', linestyle='--', label=crit_label, zorder=-10)
        crit_label = None
    plt.legend()
    plt.yscale('log')
    plt.ylim(10**(-3), 10**(3))
    plt.xscale('log')
    plt.xlabel('Height above photosphere ($R_{\odot}$)')
    plt.ylabel('Solar Wind Velocity (km/s)')




    plt.savefig(f'{this_dir}/tempest_cr{cr}.png')
    # plt.show()
    return zx, zTR, zcrit, u_miranda, u_prospero, rho


def_file = "fluxon-data/zephyr_2007_2013.sav"

def load_zephyr(file = def_file):
    from scipy.io import readsav
    data = readsav(file)
    # nz ()
    # nmods ()
    # model_year (319,)
    # tag1 (319,)
    # tag2 (319,)
    # rx (1300,)
    # rho (319, 1300)
    # uu (319, 1300)
    # valf (319, 1300)
    # t (319, 1300)
    # br (319, 1300)
    # [print(k, data[k].shape) for k in data.keys()]

    return data



def write_interpolated_file(filename, tempest_file):
    print("\n\tMaking Tempest File in Python.")
    #load the table from disk in pandas
    df = pd.read_csv(filename, delim_whitespace=True)
    df_new = df.loc[:, ['fnum', 'radius', 'b_mag']]
    # thing = df_new.groupby("fnum")

    fine_radial_distance = np.logspace(-1.25,np.log10(220), 100)
    # print(fine_radial_distance)
    fig, ax = plt.subplots(1)
    ax.set_ylabel("B_mag [Gauss]")
    ax.set_xlabel("Height above Photosphere [Rsun]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlim(10**(-2), 10**(1))
    # ax.set_ylim(10**(-2), 10**(1))
    # ax.set_ylim(4, None)

    n_model = len(df_new.groupby('fnum'))

    ax.set_title(f"{n_model} Fluxon models")

    b_mag_list = []
    first_label = "FLUX"
    doplot=True
    for group_name, group_data in df_new.groupby('fnum'):

        # Interpolate values onto the finer grid
        # interpolated_b_mag = interp1d(group_data['radius']-1, group_data['b_mag'], kind='linear', fill_value="extrapolate")(fine_radial_distance)
        interpolated_b_mag = interp1d(group_data['radius']-1, group_data['b_mag'], kind='linear', bounds_error=False)(fine_radial_distance)

        if doplot: ax.plot(fine_radial_distance, interpolated_b_mag, 'b:', label=first_label)

        b_mag_list.append(interpolated_b_mag)
        first_label = None

    if doplot:
        cr = filename.split("/data")[1].split("/")[1][-4:]
        fig.suptitle(f"FLUX vs ZEPHYR for CR {cr}")
        zephyr = load_zephyr()
        ax.plot(zephyr['rx']-1, zephyr['br'][0], 'r:', label='Zephyr')
        ax.plot(zephyr['rx']-1, zephyr['br'].T, 'r:')
        plt.legend()
        # plt.show(block=True)
        outfile = filename.split("/data")[0] + f"/imgs/mag/cr{cr}.png"
        print(f"Mag Plotting {outfile}")
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))
        plt.savefig(outfile)
        plt.close(fig)


    with open(tempest_file, 'w') as f:
        first_line = "\t" + str(len(b_mag_list)) + " " + str(len(fine_radial_distance))
        f.write(first_line)
        for rr in fine_radial_distance:
            f.write(f"\n{rr:0.8e}")

        for i in range(len(b_mag_list)):
            f.write(f"\n{i+1}  {len(fine_radial_distance)}")
            for bb in b_mag_list[i]:
                f.write(f"\n{bb:0.8e}")




def reinterpolate_velocity(tempest_file, original_data_file):
    # Load the common grid and interpolated data

    fine_radial_distance, zTR, zcrit, interpolated_velocity_miranda, interpolated_velocity_prospero, interpolated_rho = load_tempest(tempest_file)
    interpolated_velocity = interpolated_velocity_prospero

    # Load the original data for radius reference
    df_original = pd.read_csv(original_data_file, delim_whitespace=True)
    grouped = df_original.groupby('fnum')

    # Prepare DataFrame to store the results
    results = pd.DataFrame()

    # Reinterpolate for each fnum
    for fnum, group in grouped:
        original_radius = group['radius'].values
        interpolator = interp1d(fine_radial_distance, interpolated_velocity[fnum-1], kind='linear', fill_value='extrapolate')
        reinterpolated_velocity = interpolator(original_radius)

        temp_df = pd.DataFrame({
            'fnum': fnum,
            'radius': original_radius,
            'velocity': reinterpolated_velocity
        })
        results = pd.concat([results, temp_df], ignore_index=True)

    # Optionally, return or save the DataFrame
    results.to_csv(original_data_file.replace('bmag_all.dat', 'wind_tempest_reinterpolated.dat'), sep=' ', index=False)
    return results

def plot_reinterpolated_velocity(original_data_file, ax=None, show=False):

    data_file = original_data_file.replace('bmag_all.dat', 'wind_tempest_reinterpolated.dat')

    # Load the data
    df = pd.read_csv(data_file, delim_whitespace=True)

    # Group by 'fnum' to handle each line separately
    grouped = df.groupby('fnum')

    # Create a plot
    if ax is None:
        fig, ax = plt.subplots(1)

    # Loop over each group and plot
    first_label = "TEMPEST"
    for name, group in grouped:
        plt.plot(group['radius'], group['velocity'], "m", label=first_label)
        first_label = None

    # Adding title and labels
    # plt.title('Reinterpolated Velocity Profiles')
    plt.xlabel('Radius')
    plt.ylabel('Velocity')

    # Adding legend
    plt.legend()

    if show:
        # Show the plot
        plt.show()



from fluxpipe.helpers.pipe_helper import configurations

def parse_args():
    # Create the argument parser
    # print("\n\tMaking Tempest File in Python.")
    configs = configurations()
    import argparse
    parser = argparse.ArgumentParser(description=
                            'This script plots the expansion factor of the given radial_bmag_all.dat')
    parser.add_argument('--cr',     type=int, default=configs['rotations'][0],    help='Carrington Rotation')
    parser.add_argument('--dat_dir',type=str, default=configs["data_dir"],        help='data directory')
    parser.add_argument('--batch',  type=str, default=configs["batch_name"],      help='select the batch name')
    parser.add_argument('--nwant',  type=int, default=configs["fluxon_count"][0], help='magnetogram file')
    parser.add_argument('--show',   type=int, default=0)
    parser.add_argument('--file',   type=str, default=None)
    parser.add_argument('--adapt',  type=int, default=configs["adapt"],           help='Use ADAPT magnetograms')
    args = parser.parse_args()
    filename = f'{args.dat_dir}/batches/{args.batch}/data/cr{args.cr}/wind/cr{args.cr}_f{args.nwant}_radial_bmag_all.dat'
    print(filename.split("/data")[1].split("/")[1][-4:])
    configs = configurations(args=args)
    return filename, configs


def run_tempest():
    filename, configs = parse_args()

    # print(filename)
    tempest_file = filename.replace("all.dat", "tempest.dat")
    if not os.path.exists(tempest_file):
        print(f"Tempest file does not exist. Creating one at {tempest_file}")
        write_interpolated_file(filename, tempest_file)
        main(tempest_file)
        reinterpolate_velocity(tempest_file, filename)
    else:
        print(f"Tempest file exists at {tempest_file}")
    plot_tempest(tempest_file)

if __name__ == '__main__':
    run_tempest()
