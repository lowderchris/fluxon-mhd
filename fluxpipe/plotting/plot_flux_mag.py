from fluxpipe.science.read_world import rdworld
from fluxpipe.helpers.pipe_helper import configurations
configs = configurations()

# Path: fluxon-mhd/fluxpipe/fluxpipe/plotting/plot_flux_mag.py

world_path = "fluxon-data/batches/tempest7/data/cr2150/world/cr2150_f1000_hmi.flux"
world = rdworld(world_path)


import matplotlib.pyplot as plt
import numpy as np

print((world.fx.fl))
# print(len(world.fx))


plt.plot(np.asarray(world.fx.fc0), 'o', label='fc0')
plt.plot(np.asarray(world.fx.fc1), 'o', label='fc1')
plt.plot(np.asarray(world.fx.fl), 'o', label='fl')

plt.legend()
plt.show()
# a=1