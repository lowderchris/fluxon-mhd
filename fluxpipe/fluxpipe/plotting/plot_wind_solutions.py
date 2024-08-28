import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from fluxpipe.helpers.pipe_helper import configurations
from scipy.interpolate import interp1d
import re
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

def plot_full_velocity_profiles(directory=None):
    print("Plotting full velocity profiles...")
    pattern = "results_*_full_velocity.dat"
    configs = configurations()

    # Extract configuration settings
    batch = configs.get("batch_name")
    dat_dir = configs.get("data_dir")
    directory = directory or f'{dat_dir}/batches/{batch}/data/'
    # list all of the directories in the data directory
    rotations = os.listdir(directory)
    # filter out the directories that are not of the form cr#
    rotations = [d[-4:] for d in rotations if re.match(r'cr\d+', d)]

    # Use regular expression to find the integer value following 'cr'
    # match = re.search(r'/cr(\d+)/', directory)
    # cr_value = int(match.group(1))
    file = 'cr{}/wind/full_velocity_profiles/'
    rot_dirs = [directory + file.format(cr) for cr in rotations]

    for rot, cr_value in zip(rot_dirs, rotations):
        file_pattern = os.path.join(rot, pattern)
        filenames = glob.glob(file_pattern)

        if not filenames:
            print(f"No files found in {directory} matching pattern {pattern}")
            continue

        plot_lines = False
        for plot_lines in [True, False]:

            if True:
                fig, (ax, ax2) = plt.subplots(1,2, sharey='all', figsize=(10, 6),gridspec_kw={"width_ratios":[3,1]})
            else:
                fig, ax = plt.subplots()
                fig, ax2 = plt.subplots()

            cut_loc = 21.5
            hist_box = dict()
            methods = []

            # Define a common grid for interpolation
            common_grid = np.logspace(-2, 2.3, 100)  # Example grid from 10^-2 to 10^2 with 100 points

            for id, filename in enumerate(filenames):
                color = f"C{id}"
                df = pd.read_csv(filename)
                method_name = os.path.basename(filename).split('_')[1]  # Extract method name from filename
                hist_box[method_name] = []
                methods.append(method_name)

                interpolated_data = []

                fluxons = df.groupby('fluxon_position')

                for fluxon_id, group in fluxons:
                    group = group.sort_values(by='radius')  # Ensure the group is sorted by radius

                    cut_value = interp1d(group['radius'] - 1, group['velocity'], fill_value="extrapolate")(cut_loc)
                    hist_box[method_name].append(cut_value)

                    interp_vel = interp1d(group['radius'] - 1, group['velocity'], bounds_error=False, fill_value=np.nan)(common_grid)
                    interpolated_data.append(interp_vel)

                interpolated_data = np.array(interpolated_data)

                if plot_lines:
                    for data in interpolated_data:
                        ax.plot(common_grid, data, color=color, alpha=0.75,
                                zorder=-10 if method_name == "wsa" else 100 if method_name == "tempest" else None)
                else:
                    low_vel = np.nanpercentile(interpolated_data, 99., axis=0)
                    high_vel = np.nanpercentile(interpolated_data,1., axis=0)

                    # ax.fill_between(common_grid, low_vel, high_vel, color=color, alpha=0.1, zorder=-10 if method_name == "wsa" else 100 if method_name == "tempest" else None)

                    # Calculate mean and std of velocity as a function of radius
                    mean_velocity = np.nanmean(interpolated_data, axis=0)
                    std_velocity = np.nanstd(interpolated_data, axis=0)

                    # Plot mean and std
                    ax.errorbar(common_grid, mean_velocity, yerr=std_velocity, fmt='-', color=color, label=method_name, alpha=0.7, lw=4,
                                zorder=-10 if method_name == "wsa" else 100 if method_name == "tempest" else None)

            zephyr = load_zephyr()

            # Interpolating and plotting Zephyr data on the common grid
            ztop = np.nanmax(zephyr['uu'].T / 10**5, axis=1)
            zbot = np.nanmin(zephyr['uu'].T / 10**5, axis=1)
            zmean = np.nanmean(zephyr['uu'].T / 10**5, axis=1)
            zstd = np.nanstd(zephyr['uu'].T / 10**5, axis=1)

            # Interpolate Zephyr data on common grid
            zephyr_interp_top = interp1d(zephyr['rx'] - 1, ztop, bounds_error=False, fill_value=np.nan)(common_grid)
            zephyr_interp_bot = interp1d(zephyr['rx'] - 1, zbot, bounds_error=False, fill_value=np.nan)(common_grid)
            zephyr_interp_mean = interp1d(zephyr['rx'] - 1, zmean, bounds_error=False, fill_value=np.nan)(common_grid)
            zephyr_interp_std = interp1d(zephyr['rx'] - 1, zstd, bounds_error=False, fill_value=np.nan)(common_grid)

            # ax.fill_between(common_grid, zephyr_interp_bot, zephyr_interp_top, color='purple', alpha=0.1)
            ax.errorbar(common_grid, zephyr_interp_mean, yerr=zephyr_interp_std, fmt='-', color='purple', label='Zephyr_2007/13', alpha=0.7, zorder=-10, lw=5)

            # Collect interpolated values for Zephyr
            zephyr_method_name = 'Zephyr_2007_2013'
            methods.append(zephyr_method_name)
            hist_box[zephyr_method_name] = []
            for line in zephyr['uu']:
                cut_value = interp1d(zephyr['rx'] - 1, line / 10**5, fill_value="extrapolate")(cut_loc)
                hist_box[zephyr_method_name].append(cut_value)

            bins = np.linspace(0, 10**3, 20)

            for num, key in enumerate(methods):
                samples = hist_box[key]
                handle = ax2.hist(samples, bins=bins, alpha=0.5, label=key, color='purple' if '2007' in key else f"C{num}",
                                histtype='stepfilled', density=True, orientation='horizontal', align='mid')
                handle = ax2.hist(samples, bins=bins, alpha=0.95, label=None, color='k' if '2007' in key else f"C{num}",
                                histtype='step', lw=4, density=True, orientation='horizontal', align='mid')
            ax2.legend(fontsize='small')

            ax.set_xlabel('z = r/$R_\odot$ - 1')
            ax.set_ylabel('Velocity [km/s]')
            ax.set_xscale('log')
            ax.set_ylim(0, 10**3)
            ax.set_xlim(1e-2, 220)
            ax.set_title(f'Full Velocity Profiles Comparison CR {cr_value}')
            ax2.set_title(f'Hists at z={cut_loc}')
            ax.axvline(cut_loc, color='k', ls='--', alpha=0.95, lw=3)

            plt.tight_layout()

            lab = "lines_" if plot_lines else "filled_"
            outfile = f'{dat_dir}/batches/{batch}/imgs/wind/full_velocity_profiles/{lab}velocity_profile_comparison_cr{cr_value}.png'

            if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))

            plt.savefig(outfile)
            print(f"Saved plot to {outfile}")
            plt.close(fig)

if __name__ == "__main__":
    # import sys

    # if len(sys.argv) != 2:
    #     print("Usage: python plot_full_velocity_profiles.py <directory>")
    #     sys.exit(1)

    # directory = None or sys.argv[1]
    plot_full_velocity_profiles()
