import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import numpy as np


def data_shift(points):
    x = np.deg2rad(np.mod(points[0]+180, 360))
    y = np.sin(np.deg2rad(-points[1]+180))
    return x, y

# Read data from CSV
data = pd.read_csv('data.csv', header=None)

# Extract the image path
image_path = data.iloc[0, 4].split('=')[1]

# Load the FITS file
with fits.open(image_path) as hdul:
    img_data = hdul[0].data

# Normalize the image data for display purposes
normalized_img = (img_data - np.min(img_data)) / (np.max(img_data) - np.min(img_data))

# Drop the first row which contains the image path
data = data.drop(0)

# Convert coordinates to numeric
data[[0, 1, 2, 3]] = data[[0, 1, 2, 3]].apply(pd.to_numeric)

# Extract original and discovered points
original_points = data[[0, 1]].values
discovered_points = data[[2, 3]].values

# Plotting
plt.imshow(normalized_img, cmap='viridis', origin='lower', extent=[0, 2*np.pi, -1, 1])  # Adjust color map if needed

# Draw arrows from original to discovered points
for original, discovered in zip(original_points, discovered_points):
    orig = data_shift(original)
    dis = data_shift(discovered)
    dff = np.array(dis) - np.array(orig)
    plt.scatter(*orig, color="red")
    plt.scatter(*dis, color="blue")
    plt.arrow(*orig, *dff, color='green', head_width=0.025, length_includes_head=True)

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Footpoints and Nearest Boundary')
plt.show()
