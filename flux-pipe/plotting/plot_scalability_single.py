import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, StrMethodFormatter
import pandas as pd
# import statsmodels.api as sm

run_name = "fluxon r3"
datdir = "/Users/cgilbert/vscode/fluxon-data/"
file_path = os.path.join(datdir, f"{run_name}/pfss_time.txt")


## Helper functions to parse the output file
def convert_value(value):
    """Convert a string to an int or float if possible, otherwise return the string."""
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value.strip()

def parse_line(line):
    """Parse a line of the output file into a dictionary."""
    key_values = line.strip().split(',')
    # print(key_values)
    parsed_dict = {}
    for key_value in key_values:
        if ':' in key_value:
            key, value = key_value.split(':')
            parsed_dict[key.strip()] = convert_value(value)
    return parsed_dict

def load_data(file_path):
    """Load the data from the file into a pandas DataFrame."""
    data = []
    print("\n", file_path, '\n')
    with open(file_path, 'r') as file:
        for line in file.readlines():
            data.append(parse_line(line))
    return pd.DataFrame(data)



## Real Code

# Load the Array
df = load_data(file_path)
df = df.drop(df.index[0])
print(df)

cr = df['cr'].unique()[0]

n_reductions = len(df["r"].unique())
ny_plots = min(2, n_reductions)
nx_plots = min(int(np.ceil(n_reductions/ny_plots)), n_reductions)
fig, axarray = plt.subplots(nx_plots, ny_plots, figsize=(12, 8), sharex="all", sharey="all")
if type(axarray) != np.ndarray:
    axarray = np.array(axarray)
flaxarray = axarray.flatten()


if len(df["r"].unique()) == 1:
    reduce = int(df["r"].unique()[0])
else:
    reduce = 0

to_plot = ["n_want", "TrOpen", "n_out", 'ttime', "Success"]
independent = "n_actual"
for (ax, reduction, rx, ry) in zip(flaxarray, df["r"].unique(), df['rx'].unique(), df['ry'].unique()):
    ax.set_title(f"CR{int(cr)}, Resolution: ({int(rx)}, {int(ry)})")
    ax.set_yscale('log')
    ax.set_xscale('log')

    # format = ScalarFormatter()
    # format.set_scientific(False)
    # format.ticklabel_format(style='plain')
    # ax.xaxis.set_major_formatter(format)
    # plt
    # ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    # ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    # ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=False))


    ax.set_xlabel("Actual Number of Fluxons")
    ax.set_ylabel("Values", labelpad=-5)
    # ax.yaxis.label.set_offset(10)
    for ii, line_name in enumerate(to_plot):
        this_line = df[df["r"] == reduction][line_name]
        abscissa  = df[df["r"] == reduction][independent]

        # Filter NaN values
        valid_indices = np.logical_and(~np.isnan(abscissa), ~np.isnan(this_line))
        abscissa_valid = abscissa[valid_indices]
        this_line_valid = this_line[valid_indices]
        if abscissa_valid.any():
            coefficients = np.polyfit(abscissa_valid, this_line_valid, 1)
            fit_line = np.poly1d(coefficients)(abscissa)
            round_coeff = [round(ccc, 2) for ccc in coefficients]
            this_label = f"{line_name}: {round_coeff[0]}"
            color = f"C{ii}"

            if line_name in ["Success"]: 
                this_line_valid[this_line_valid <= 0.1] = 0.1
                # ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label)
                # the same plot but stair-stepped
                ax.step(abscissa_valid, this_line_valid, color=color, where='mid', marker='o', label="Success")
                # find the index of the first non-"1" value in this_line_valid
                non_one_index = next((i for i, val in enumerate(this_line_valid) if val != 1), None)
                non_one_loc = abscissa_valid[non_one_index]
                ax.axvline(non_one_loc, color="k", ls=":")

                
            elif line_name in ["ttime"]: 
                ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label, marker='o')
            else:
                ax.plot(abscissa_valid, this_line_valid, color=color, label=this_label, marker='o')
                this_line_valid[this_line_valid <= 0.1] = 0.1
                fit_line[fit_line <= 0.1] = 0.1
                # fit_line[fit_line == 0] = 0.1
                ax.plot(abscissa, fit_line,  color=color, ls=":")
    ax.axhline(100, color="k", ls="--", label="100 Fluxons")
    # ax.legend(loc="center right") 
    ax.legend(loc=(0.6, 0.25 ))
# set the ylim to go from 1900 to the current ymax
# ax.set_xlim(10**2, ax.get_xlim()[1]) # this is the magic

ax.set_ylim((5*10**(-2), 2.1*(10**4)))
ax.set_xlim((80, 15000))

# line that changes the ylim to be from 1 to the current ymax
fig.set_size_inches(5, 5.5)
# fig.suptitle("Scalability of Fluxon Algorithm")
# Fit line
plt.tight_layout()

img_path = os.path.join(datdir, f"{run_name}/imgs/scalability_{reduce}.png")
img_path_pdf = os.path.join(datdir, f"{run_name}/imgs/scalability_{reduce}.pdf")
plt.tight_layout()

plt.savefig(img_path)
plt.savefig(img_path_pdf)
# plt.show()
plt.close(fig)




# regression.summary
# data = pd.DataFrame(index=x, data={'y': y, 'trend': trend})
# df.plot(x="n_actual", y=["n_pfss", "n_out", "TrOpen", "ttime"], ax=ax) # add kwargs for title and other layout/design aspects

# def load_data(file_path, dependent, independent):
#     data = []
#     print("\n", file_path, '\n')
#     with open(file_path, 'r') as file:
#         for line in file.readlines():
#             data.append(parse_line(line))
#     df = pd.DataFrame(data)
#     X = sm.add_constant(df[independent])
#     model = sm.OLS(df[dependent], X)
#     regression = model.fit()
#     trend = regression.predict(X)
#     return pd.DataFrame({'data': data, 'trend': trend})


# Example usage
# file_path = 'data.txt'  # Replace with your actual file path




# def load_data(file_path):
#     data = []
#     print("\n", file_path, '\n')
#     with open(file_path, 'r') as file:
#         for line in file.readlines():
#             data.append(parse_line(line))
#     regression = pd.ols(y=dependent, x=independent)
#     trend = regression.predict(beta=regression.beta, x=df[dependent]) 
#     return pd.DataFrame(data ,'trend': trend)

# # import pandas as pd

# # Define the file path
# file_path = f'{datdir}/fluxon/pfss_time copy 2.txt'

# import csv

# # Initialize empty lists for each column
# rez = []
# x,y = [],[]
# time = []
# time_kPix = []
# trOpen = []
# trClosed = []
# fail = []
# fluxons = []
# steps = []

# # Read the regularized table from a CSV file
# with open(file_path, 'r') as file:
#     reader = csv.reader(file, delimiter='|')
#     next(reader)  # Skip header row

#     for row in reader:
#         row = [item.strip() for item in row]  # Remove leading/trailing spaces
#         # Extract values from each column
#         r, rez_value, time_value, time_kPix_value, trOpen_value, trClosed_value, \
#         fail_value, fluxons_value, steps_value = row

#         # Append values to corresponding lists
#         rez.append(rez_value)
#         time.append(float(time_value))
#         time_kPix.append(float(time_kPix_value))
#         trOpen.append(int(trOpen_value))
#         trClosed.append(int(trClosed_value))
#         fail.append(int(fail_value))
#         fluxons.append(int(fluxons_value))
#         steps.append(int(steps_value))

# # Print the loaded variables

# x,y = [(x) for x,y in rez], [(y) for x,y in rez]
# print(x,y)
# print(rez)
# print(time)
# print(time_kPix)
# print(trOpen)
# print(trClosed)
# print(fail)
# print(fluxons)
# print(steps)





# # # Define the column names
# # column_names = ['r', 'rez', 'time', 'TrOpen', 'TrClosed', 'Fail', 'Fluxons', 'Steps', 'Stiff']

# # # Initialize an empty list to store the data
# # data = []

# # # Read the file line by line
# # with open(file_path, 'r') as file:
# #     for line in file:
# #         # Split the line into individual values
# #         values = line.strip().split(',')
        
# #         # Convert the values to the appropriate data types
# #         r = int(values[0].split('=')[1].strip())
# #         rez = tuple(map(int, values[1].split('=')[1].strip()[1:-1].split(',')))
# #         time = float(values[2].split(':')[0].strip().split(' ')[0])
# #         TrOpen = int(values[3].split(':')[1].strip())
# #         TrClosed = int(values[4].split(':')[1].strip())
# #         Fail = int(values[5].split(':')[1].strip())
# #         Fluxons = int(values[6].split(':')[1].strip().split(' ')[0])
# #         Steps = int(values[6].split(':')[2].strip().split(' ')[0])
# #         Stiff = values[7].strip().split(' ')[-1] if len(values) > 7 else None
        
# #         # Append the values to the data list
# #         data.append([r, rez, time, TrOpen, TrClosed, Fail, Fluxons, Steps, Stiff])

# # # Create a pandas DataFrame from the data
# # df = pd.DataFrame(data, columns=column_names)
















# # # read_fr.py
# # # Script to parse fluxon output file of the field expansion factor


# # # file output format: fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1
# # # fluxon id, beginning coords, end coords, beginning mag, end mag, area maybe?

# # import numpy as np
# # import matplotlib as mpl
# # mpl.use("qt5agg")
# # import matplotlib.pyplot as plt
# # import argparse
# # from astropy.io import fits
# # from magnetoget import load_fits_magnetogram

# # print("Plotting Bmag...", end="")
# # # create the argument parser
# # parser = argparse.ArgumentParser(description='This script plots the expansion factor of the given radial_fr.dat')
# # parser.add_argument('--cr', type=int, default=0, help='Carrington Rotation')
# # parser.add_argument('--dat_dir', type=str, default='/Users/cgilbert/vscode/fluxon-data', help='data directory')
# # parser.add_argument('--show', type=int, default=0)
# # args = parser.parse_args()
# # from magnetoget import load_magnetogram_params
# # (hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(args.dat_dir)
# # CR = args.cr or cr or 2183

# # filename = f'{args.dat_dir}/fluxon/cr{CR}/wind/radial_bmag.dat'

# # # Load the dat file
# # arr = np.loadtxt(filename).T
# # fid, phi0, theta0, phi1, theta1, br0, br1, ar0, ar1 = arr
# # nfluxon = arr.shape[1]


# # # Convert coords to correct coords
# # ph0, th0 = phi0+np.pi, -(theta0-(np.pi/2))
# # ph1, th1 = phi1+np.pi, -(theta1-(np.pi/2))


# # # Do some data manipulation
# # br0_max = np.nanmax(br0) or 0.25
# # br1_max = np.nanmax(br1) or 0.25
# # ar0_max = np.nanmax(ar0) or 0.25
# # ar1_max = np.nanmax(ar1) or 0.25

# # skew = 5**2
# # power = 1
# # b0 = skew*(6*br0/br0_max)**power
# # b1 = skew*(4*br1/br1_max)**power
# # a0 = skew*(6*ar0/ar0_max)**power
# # a1 = skew*(4*ar1/ar1_max)**power


# # # Plot the Data
# # fig, (ax0, ax1) = plt.subplots(2)

# # magnet = load_fits_magnetogram()
# # ax0.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')
# # ax1.imshow(magnet, cmap='gray', interpolation=None, origin="lower", extent=(0,2*np.pi,-1,1), aspect='auto')

# # sc00 = ax0.scatter(ph0, th0, c=br0, s = b0, cmap="winter", alpha=0.75, label='B(1.0R$_\odot$)')
# # sc01 = ax0.scatter(ph1, th1, c=br1, s = b1, cmap="autumn", alpha=0.75, label='B(21.5R$_\odot$)', marker='s')
# # sc10 = ax1.scatter(ph0, th0, c=ar0, s = a0, cmap="winter", alpha=0.75, label='A(1.0R$_\odot$)')
# # sc11 = ax1.scatter(ph1, th1, c=ar1, s = a1, cmap="autumn", alpha=0.75, label='A(21.5R$_\odot$)', marker='s')

# # cbar01 = fig.colorbar(sc01, ax=ax0)
# # cbar00 = fig.colorbar(sc00, ax=ax0)
# # cbar11 = fig.colorbar(sc11, ax=ax1)
# # cbar10 = fig.colorbar(sc10, ax=ax1)

# # cbar00.set_label("B(1.0R$_\odot$)  [Gauss] ")
# # cbar01.set_label("B(21.5R$_\odot$) [Gauss]")
# # cbar10.set_label("A(1.0R$_\odot$)  [m$^2$]")
# # cbar11.set_label("A(21.5R$_\odot$) [m$^2$]")

# # ax0.set_title(F"Fluxon Magnetic Field Strength for CR {args.cr}")
# # ax1.set_title(F"Fluxon Area for CR {args.cr}")

# # for ax in (ax0, ax1):
# #     ax.set_xlabel('Longitude (Radians)')
# #     ax.set_ylabel('Sine latitude')
# #     ax.set_ylim((-1.5,1.1))
# #     ax.axhline(-1, c='lightgrey', zorder=-10)
# #     ax.axhline( 1, c='lightgrey', zorder=-10)
# #     ax.axvline(0, c='lightgrey', zorder=-10)
# #     ax.axvline(2*np.pi, c='lightgrey', zorder=-10)
# #     ax.legend(loc="lower right")

# # fig.set_size_inches((8,8))
# # plt.tight_layout()
# # pngname = filename.replace(".dat", f"_{len(ph0)}_{len(ph1)}.png")
# # # pngname = filename.replace(".dat", ".png")
# # plt.savefig(pngname)
# # if args.show:
# #     plt.show()
# # plt.close(fig)
# # print("Done!")





# # # ax0.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*br0/br0_max)**3, alpha=0.75, label='B(1.0R$_\odot$)')
# # # ax0.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*br1/br1_max)**3, alpha=0.75, label='B(21.5R$_\odot$)')
# # # ax1.scatter(phi0%(2*np.pi), np.sin(theta0), s=(1+5*ar0/ar0_max)**3, alpha=0.75, label='A(1.0R$_\odot$)', marker='s')
# # # ax1.scatter(phi1%(2*np.pi), np.sin(theta1), s=(1+5*ar1/ar1_max)**3, alpha=0.75, label='A(21.5R$_\odot$)', marker='s')