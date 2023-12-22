import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import numpy as np

# Read data from CSV
data = pd.read_csv('data.csv')

# Extract the image path
image_path = data['image_path'].iloc[0]

# Load the FITS file
with fits.open(image_path) as hdul:
    img_data = hdul[0].data

# Normalize the image data for display purposes
normalized_img = (img_data - np.min(img_data)) / (np.max(img_data) - np.min(img_data))

# Filter original and discovered points
original = data[data['type'] == 'original']
discovered = data[data['type'] == 'discovered']

# Plotting
plt.imshow(normalized_img, cmap='gray', origin='lower')  # Adjust color map if needed
plt.scatter(original['x'], original['y'], color='red', label='Original Point')
plt.scatter(discovered['x'], discovered['y'], color='blue', label='Discovered Point')

plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Original and Discovered Points')
plt.legend()
plt.show()
