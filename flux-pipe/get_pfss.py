"""
Fluxon PFSS mapping
=======================
Generates a fluxon mapping from input GONG-sourced pfss coronal field solution
"""


datdir = "/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/"


###############################################################################
# First, import required modules
# import matplotlib as mpl
# mpl.use("qt5agg")
import os
import sys
import numpy as np
from magnetoget import load_magnetogram_params

# Check for command line arguments
if len(sys.argv)>1: 
    force_plot = force_trace = True if int(sys.argv[1])>0 else False
else:
    force_plot = force_trace = False

try:
    batch = sys.argv[2]
except IndexError:
    batch = "fluxon_paperfigs_2"

###############################################################################
# Load the parameters, basically just which magnetogram is used
print("\n -->Loading Parameters...")
(hdr, cr, fname, adapt, doplot, reduce) = load_magnetogram_params(datdir)
elapsed = 0



###############################################################################
# Load the fits file and format the data and header
from pfss_funcs import load_and_condition_fits_file
br_safe, fits_path = load_and_condition_fits_file(fname, datdir)



###############################################################################
# Do the PFSS mapping
from pfss_funcs import load_pfss, compute_pfss, fits_path_to_pickle_path
pickle_dir = os.path.join(datdir, "pfss")
if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
pickle_path = os.path.join(pickle_dir, f"pfss_cr{cr}_r{reduce}.pkl")

output = load_pfss(pickle_path)
if not output:
    output, elapsed = compute_pfss(br_safe, pickle_path) #, nrho, rss)



###############################################################################
# Plot the fluxon locations
from pfss_funcs import plot_fluxon_locations
f_lon, f_lat, f_sgn = plot_fluxon_locations(br_safe, cr, datdir, fits_path, reduce, force_plot, batch)
n_flux = len(f_sgn)



###############################################################################
# Trace pfss field lines
from pfss_funcs import trace_lines
skip_num = 'x'
timeout_num = 'x'
open_path = datdir + f'{batch}/cr' + cr + '/floc/' + 'floc_open_cr'+cr+'.dat'
closed_path = datdir + f'{batch}/cr' + cr + '/floc/' + 'floc_closed_cr'+cr+'.dat'

print("\n -->Tracing Open and Closed Fluxons...")

if not os.path.exists(open_path) or force_trace:
    fl_open, fl_closed, skip_num, timeout_num, flnum_open, flnum_closed = trace_lines(output, f_lon, f_lat, f_sgn, open_path, closed_path)
else:
    fl_open = np.loadtxt(open_path)
    fl_closed = np.loadtxt(closed_path)
    flnum_open = len(np.unique(fl_open[:, 0]))+1
    flnum_closed = 2*len(np.unique(fl_closed[:, 0]))
    print("\tSkipped! floc dat files already exist.")
    print(f"\tFootpoints:\t Open: {flnum_open}, Closed: {flnum_closed}, Total: {flnum_open+flnum_closed}")
    # print(f"\tFluxons:\t Open: {flnum_open}, Closed: {flnum_closed//2}, Total: {flnum_open+flnum_closed//2}")





###############################################################################
# Record stats in the output file
shp = br_safe.data.shape
pix = shp[0]*shp[1]
timefile = datdir + f'{batch}/pfss_time.txt'
with open(timefile, 'a+') as f:
    # a good name for the variable
    elap =f"\ncr: {cr}, r: {reduce}, rx: {shp[0]}, ry: {shp[1]}, pf_elp: {elapsed:0>3.3f}, t_kpix: {1000*elapsed/pix:0.3f}"
    nlines = f"TrOpen: {flnum_open}, TrClosed: {flnum_closed}, TrGood: {flnum_open+flnum_closed}, TrFail: {skip_num+timeout_num}, "
    f.write(f"{elap}, {nlines}")
    # f.write(")




















# import pandas as pd

# # assuming shp and reduce are already defined

# # create a dictionary of the data to be appended
# data = {
#     'r': [reduce],
#     'rez': [shp],
#     'TrOpen': [len(fl_open)],
#     'TrClosed': [len(fl_closed)],
#     'Fail': [skip_num+timeout_num],
#     'elapsed': [elapsed],
#     'pix_time': [1000*elapsed/pix],
# }

# # create a DataFrame from the data
# df = pd.DataFrame(data)

# # define the file path to append the data
# timefile = datdir + 'fluxon/pfss_time.txt'

# # append the data to the file
# df.to_csv(timefile, mode='a', header=False, index=False, sep='\t')

# # print("\n-->>>>>>>>>>>>>>\n-->>>>Main Program Complete!<<<<\n<<<<<<<<<<<<<<\n")
# print(" -->Plotting: ", bool(doplot), end="\n\n")
# def set_axes_lims(ax):
#     ax.set_xlim(0, 360)
#     ax.set_ylim(0, 180)

# # If you want to plot...
# doplot = False
# if doplot and False:
#     print("\n -->Plotting...")

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     r = 1.01
#     for i in np.arange(0,len(f_lon)):
#         x0 = SkyCoord(f_lon[i] * u.rad, f_lat[i] * u.rad, r)

#         # x0 = coords.sph2cart(r, f_lat[i], f_lon[i])
#         fl = output.trace(tracer, x0)
#         fl = fl.field_lines[0]
#         color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(fl.polarity)
#         ax.plot(fl.x / const.R_sun,
#                 fl.y / const.R_sun,
#                 fl.z / const.R_sun,
#                 color=color, linewidth=1)
#         if fl.is_open:
#             fl = fl.spherical
#             x0 = np.array(coords.sph2cart(2.45, -1*(-1*fl.lat[-1].value-90)*np.pi/180, fl.lon[-1].value*np.pi/180))
#             fl = output.trace(x0)
#             color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(fl.polarity)
#             ax.plot(fl.x / const.R_sun,
#                     fl.y / const.R_sun,
#                     fl.z / const.R_sun,
#                     color=color, linewidth=1)
#     #ax.set_title('CR'+cr+' PFSS field')
#     plt.show()

# from matplotlib.pyplot import scatter, plot, imshow, xlim, ylim, xlabel, ylabel
# from numpy import where

# if doplot:
#     scatter(fl_open[:,2], fl_open[:,1])
#     scatter(fl_open[:,5], fl_open[:,4])
#     for i in np.arange(0, len(fl_open[:,0])):
#        plot([fl_open[i,2], fl_open[i,5]],[fl_open[i,1], fl_open[i,4]],'k')

#     scatter(fl_closed[:,2], fl_closed[:,1])
#     scatter(fl_closed[:,5], fl_closed[:,4])
#     for i in np.arange(0, len(fl_closed[:,0])):
#        plot([fl_closed[i,2], fl_closed[i,5]],[fl_closed[i,1], fl_closed[i,4]],'k')

#     wp = where(f_sgn>0)[0]
#     wn = where(f_sgn<0)[0]
#     scatter(f_lon[wn], f_lat[wn])
#     scatter(f_lon[wp], f_lat[wp])

#     wp = where(fl_open[:,6]>0)[0]
#     wn = where(fl_open[:,6]<0)[0]
#     scatter(fl_open[wn,2], fl_open[wn,1])
#     scatter(fl_open[wp,2], fl_open[wp,1])

#     imshow(br, vmin=-25, vmax=25, extent=[0,360,-90,90], aspect='auto')
#     scatter(f_lon[where(f_sgn>0)]*180/np.pi, -1*f_lat[where(f_sgn>0)]*180/np.pi+90, c='white')
#     scatter(f_lon[where(f_sgn<0)]*180/np.pi, -1*f_lat[where(f_sgn<0)]*180/np.pi+90, c='black')
#     xlim([0,360])
#     ylim([-90,90])
#     xlabel('Carrington longitude')
#     ylabel('Latitude')

#     cmap = plt.get_cmap('tab10')
#     imshow(br, vmin=-20, vmax=20, extent=[0,360,-1,1], aspect='auto')
#     scatter(fl_closed[:,2], np.sin(fl_closed[:,1]*np.pi/180), s=5.0, color=cmap(0))
#     scatter(fl_closed[:,5], np.sin(fl_closed[:,4]*np.pi/180), s=5.0, color=cmap(1))
#     scatter(fl_open[:,2], np.sin(fl_open[:,1]*np.pi/180), s=5.0, color=cmap(3))
#     xlim([0,360])
#     ylim([-1,1])
#     xlabel('Carrington longitude (CR2193)')
#     ylabel('Sine latitude')

#     plt.show()


# if __name__ == "__main__":
#     sys.exec("python3 mag_runner.py")



    # An attempt to use the performance cores
    # from os import setpriority, getpid, getpriority
    # PRIO_DARWIN_THREAD  = 0b0011
    # PRIO_DARWIN_PROCESS = 0b0100
    # PRIO_DARWIN_BG      = 0x1000
    # setpriority(PRIO_DARWIN_PROCESS, 0, 0)

    # # Actually do the math
    # import subprocess
    # import time
    # before = time.time()

    # # Serialize the input object and write it to a file
    # with open('input.pkl', 'wb') as file:
    #     pickle.dump(input, file)

    # # Create a subprocess
    # command = f"import pfsspy; import pickle; pickle.dump(pfsspy.pfss(pickle.load(open('input.pkl', 'rb'))), open('{pickle_path}', 'wb'), pickle.HIGHEST_PROTOCOL)"
    # process = subprocess.Popen(['python', '-c', command])

    # # Loop to periodically check the subprocess status
    # while True:
    #     time.sleep(2)  # Wait for 2 seconds before checking again

    #     # Check if the subprocess has completed
    #     if process.poll() is not None:
    #         # Subprocess has completed
    #         break

    #     # Subprocess is still running
    #     print('.', end="")

    # # Subprocess has completed, retrieve the output
    # output = process.communicate()[0]

    # # Process the output as desired
    # print(output)

    # import multiprocessing
    # import time

    # input_object = input  # Replace with your actual input object

    # # Define a function to run in the subprocess
    # def run_pfss(input_obj, output_queue):
    #     import pfsspy

    #     output = pfsspy.pfss(input_obj)
    #     output_queue.put(output)

    # # Create a queue for receiving the output
    # output_queue = multiprocessing.Queue()

    # # Create a subprocess
    # process = multiprocessing.Process(target=run_pfss, args=(input_object, output_queue))
    # print("Starting process...")
    # process.start()
    # print("Process started!")
    # exit()

    # # Loop to periodically check the subprocess status
    # while True:
    #     time.sleep(2)  # Wait for 2 seconds before checking again

    #     # Check if the subprocess has completed
    #     if not process.is_alive():
    #         # Subprocess has completed
    #         break

    #     # Subprocess is still running
    #     # Perform any other tasks or checks here
    #     print('.', end="")

    # # Subprocess has completed, retrieve the output
    # output = output_queue.get()



    # exit()


    # import subprocess
    # # import time

    # # Create a subprocess
    # process = subprocess.Popen(['python', '-c', 'import pfsspy; output = pfsspy.pfss(input)'])

    # # Loop to periodically check the subprocess status
    # while True:
    #     time.sleep(5)  # Wait for 2 seconds before checking again
        
    #     # Check if the subprocess has completed
    #     if process.poll() is not None:
    #         # Subprocess has completed
    #         break
        
    #     # Subprocess is still running
    #     # Perform any other tasks or checks here
    #     print(".", end="")
    # # Subprocess has completed, retrieve the output
    # output = process.communicate()[0]

    # # Process the output as desired
    # # print(output)



    # import pdb
    # pdb.set_trace()
    # ax.imshow(output.bc[0], cmap='RdBu', interpolation=None, origin="lower", zorder=-5, alpha= 0.33)
    # pfss_out = output
    # ss_br = pfss_out.source_surface_br
    # # Create the figure and axes
    # # fig2 = plt.figure()
    # ax = plt.subplot(projection=ss_br)
    # fig, ax = plt.subplots()
    # # Plot the source surface map
    # ss_br.plot()
    # # Plot the polarity inversion line
    # ax.plot_coord(pfss_out.source_surface_pils[0])
    # # Plot formatting
    # plt.colorbar()
    # ax.set_title('Source surface magnetic field')

    # plt.show()
    # output.bc
    # output.bg
    # output.bunit