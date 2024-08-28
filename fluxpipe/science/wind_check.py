import numpy as np
import matplotlib as mpl

mpl.use("qt5agg")

def wind_speed(fss, theta_b):

    # Define constants (From Schonfeld 2022)
    c1 = 2/9
    c2 = 0.8
    c3 = 2.0
    c4 = 2.0
    c5 = 3.0

    v0 = 285 #km/s
    vm = 910 - v0 #km/s

    speed = v0 + (vm / (1 + fss)**c1) * (1 - c2 * np.exp(-(theta_b / c3)**c4))**c5
    return speed


expansions = np.linspace(1.0, 20, 100)
angles = np.linspace(0, 20, 1000)

themap = np.zeros((len(expansions), len(angles)))

for ii, f in enumerate(expansions):
    for jj, t in enumerate(angles):
        # print(f, t, wind_speed(f, t))
        themap[ii,jj]=wind_speed(f, t)

import matplotlib.pyplot as plt
plt.imshow(themap, extent=[0, 20, 1.0, 20], aspect='auto', origin='lower')
plt.xlabel('Theta_b')
plt.ylabel('FSS')
plt.colorbar(label="Wind Speed")
plt.title("Wind Speed as a Function of Theta_footpoint and the Expansion Factor")
plt.show()