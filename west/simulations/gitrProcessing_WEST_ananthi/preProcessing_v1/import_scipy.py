import numpy as np
from scipy.io import savemat
import matplotlib.pyplot as plt

# Generate the `coolwarm` colormap with 256 colors
cmap = plt.get_cmap('coolwarm', 256)
coolwarm_rgb = cmap(np.linspace(0, 1, 256))[:, :3]  # Get RGB values

# Save the RGB data to a .mat file
savemat('coolwarm.mat', {'coolwarm_rgb': coolwarm_rgb})

