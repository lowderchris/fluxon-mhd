import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# Generate a random set of points
num_points = 20
points = np.random.rand(num_points, 2)

# Create a synthetic image (for dimensions)
image_width, image_height = 200, 200
image = np.zeros((image_height, image_width))

# Generate grid points based on the image shape
y_indices, x_indices = np.indices(image.shape)
grid_points = np.column_stack((y_indices.ravel(), x_indices.ravel()))

# Scale points to the image dimensions
scaled_points = np.copy(points)
scaled_points[:, 0] *= image_width
scaled_points[:, 1] *= image_height

# Build a KDTree for efficient nearest-neighbor query
tree = cKDTree(scaled_points)

# Query the nearest distance for each grid point
distances, _ = tree.query(grid_points)

# Reshape the distances to the image shape
distance_array = distances.reshape(image.shape)

# Display the image
plt.imshow(distance_array, cmap='viridis')
plt.colorbar(label='Distance to Nearest Point')
plt.title('Distance to Nearest Random Point')
plt.show()
