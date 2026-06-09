import numpy as np
import os
from scipy.interpolate import interpn
from scipy.io import loadmat
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

# Clear all open figures
plt.close('all')


# Load .mat file using scipy.io for MATLAB data
mat_data = loadmat('west_geom.mat')
planes = mat_data['planes']  # Replace 'planes' with the actual key in your .mat file
centroid = mat_data['centroid']  # Replace 'planes' with the actual key in your .mat file
area = mat_data['area']  # Replace 'planes' with the actual key in your .mat file
# centroid = mat_data['norm_vec']  # Replace 'planes' with the actual key in your .mat file
# centroid = np.loadtxt('centroid.csv', delimiter=',')
# norm_vec = np.loadtxt('norm_vec.csv', delimiter=',')

#import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Assuming `planes` and `centroid` are loaded previously in the script
nP = int(1e6)
radius = 0.06256

# Initialize arrays for material properties
materialZ = np.zeros(len(planes))
surfs = np.zeros(len(planes))

# Define HeliconRefined set based on conditions
HeliconRefined = np.where(
    (np.abs(centroid[:, 2] - 0.16) <= 2.94) &
    (centroid[:, 0] <= 0.3583) &
    (centroid[:, 0] >= -0.3583)
)[0]

materialZ[HeliconRefined] = 13
surfs[HeliconRefined] = 1

# Plotting GITR Geometry
plotSet_no_surfs = np.setdiff1d(np.arange(len(planes)), np.where(surfs)[0])

# Calculate dynamic limits based on the data range
x_min, x_max = planes[:, [0, 3, 6]].min(), planes[:, [0, 3, 6]].max()
y_min, y_max = planes[:, [1, 4, 7]].min(), planes[:, [1, 4, 7]].max()
z_min, z_max = planes[:, [2, 5, 8]].min(), planes[:, [2, 5, 8]].max()

# Increase Z-dimension limit as specified
z_max = max(z_max, 5)

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot non-surface planes
for i in plotSet_no_surfs:
    X = [planes[i, 0], planes[i, 3], planes[i, 6]]
    Y = [planes[i, 1], planes[i, 4], planes[i, 7]]
    Z = [planes[i, 2], planes[i, 5], planes[i, 8]]
    verts = [list(zip(X, Y, Z))]
    poly = Poly3DCollection(verts, color=[0.7, 0.7, 0.7], alpha=0.3, edgecolor='none')
    ax.add_collection3d(poly)

# Set axis labels and dynamic limits
ax.set_title("GITR Geometry")
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_zlim(z_min, z_max)

# Plot surface planes in blue
plotSet_surfs = np.where(surfs)[0]
for i in plotSet_surfs:
    X = [planes[i, 0], planes[i, 3], planes[i, 6]]
    Y = [planes[i, 1], planes[i, 4], planes[i, 7]]
    Z = [planes[i, 2], planes[i, 5], planes[i, 8]]
    verts = [list(zip(X, Y, Z))]
    poly = Poly3DCollection(verts, color='b', alpha=0.3, edgecolor='none')
    ax.add_collection3d(poly)

plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Assuming `planes` is already defined as an array where each row represents a plane with three vertices
# Each plane has 9 values: x1, y1, z1, x2, y2, z2, x3, y3, z3.

# Calculate `abcd` for each plane
num_planes = len(planes)
abcd = np.zeros((num_planes, 4))  # Initialize abcd array to store plane coefficients
plane_norm = np.zeros(num_planes)  # Array to store the norms of each plane's normal vector

for i in range(num_planes):
    # Extract points defining the plane
    p1 = planes[i, 0:3]
    p2 = planes[i, 3:6]
    p3 = planes[i, 6:9]
    
    # Calculate two vectors within the plane
    v1 = p2 - p1
    v2 = p3 - p1
    
    # Calculate the normal vector using the cross product
    normal_vector = np.cross(v1, v2)
    norm = np.linalg.norm(normal_vector)
    
    # Normalize the normal vector and store it in `abcd`
    abcd[i, :3] = normal_vector / norm
    abcd[i, 3] = -np.dot(abcd[i, :3], p1)  # Calculate the d component of the plane equation
    
    # Store the plane's normal magnitude
    plane_norm[i] = norm

# Now proceed with the rest of the code

# Initialize `inDir` array and set its values to 1
inDir = np.ones(num_planes)

# Prepare figure
fig = plt.figure(10)
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Surface Normals")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Calculate normals for each plane
l_normal = 0.01  # Length factor for normal
norm_vec = np.zeros((num_planes, 3))  # Array to store normal vectors

for i in range(num_planes):
    # Normal calculation
    normal = -abcd[i, :3] / plane_norm[i]
    normal = l_normal * normal

    # Calculate centroid (average of the three vertices in the plane)
    centroid[i] = [
        np.mean([planes[i, 0], planes[i, 3], planes[i, 6]]),
        np.mean([planes[i, 1], planes[i, 4], planes[i, 7]]),
        np.mean([planes[i, 2], planes[i, 5], planes[i, 8]])
    ]

    # Determine inDir based on dot product
    dot_product = np.dot(normal, centroid[i])
    if np.sign(dot_product) == 1:
        inDir[i] = -1

    # Apply inDir to adjust the direction of the normal
    normal = inDir[i] * normal
    norm_vec[i, :] = normal  # Store normal in norm_vec

    # Optional: Visualize the normal vector using quiver3
    # vec_factor = 1
    # ax.quiver(centroid[i, 0], centroid[i, 1], centroid[i, 2], 
    #           vec_factor * normal[0], vec_factor * normal[1], vec_factor * normal[2])

# Save centroid data to 'centroid.csv'
np.savetxt('centroid.csv', centroid, delimiter=',')
np.savetxt('norm_vec.csv', norm_vec, delimiter=',')
norm_vec = np.loadtxt('norm_vec.csv', delimiter=',')

plt.show()


# # Calculate surface properties
# inDir = np.ones(len(planes))
# # norm_vec_mag = np.linalg.norm(norm_vec, axis=1)
# # unorm_vec = norm_vec / norm_vec_mag[:, None]

# # Loop over planes to adjust normals and directions
# for i in range(len(planes)):
#     normal = -planes[i, 0:3] / np.linalg.norm(planes[i, 0:3])
#     dot_product = np.dot(normal, centroid[i])
#     if dot_product > 0:
#         inDir[i] = -1
#     normal *= inDir[i]

#     if abs(normal[2]) > 0.005:
#         planes[i, -1] = 0
#         planes[i, -2] = 0

# # Save updated centroid and norm_vec back to CSV
# np.savetxt('centroid.csv', centroid, delimiter=',')
# np.savetxt('norm_vec.csv', norm_vec, delimiter=',')
# norm_vec = np.loadtxt('norm_vec.csv', delimiter=',')

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Constants
ME = 9.10938356e-31
MI = 1.6737236e-27
Q = 1.60217662e-19
EPS0 = 8.854187e-12
amu = 18

# Geometry and profile data
print(">>>> Loading geometry and profile data")
mat_data = loadmat('west_geom.mat')
planes = mat_data['planes']
centroid = np.loadtxt('centroid.csv', delimiter=',')
norm_vec = np.loadtxt('norm_vec.csv', delimiter=',')

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot as plt

# Physical Constants
ME = 9.10938356e-31
MI = 1.6737236e-27
Q = 1.60217662e-19
EPS0 = 8.854187e-12
amu = 18

# Load profiles from NetCDF and reshape if necessary
with Dataset('profilesWEST.nc', 'r') as nc:
    R = nc.variables['x'][:]  # Radial grid points
    z = nc.variables['z'][:]  # Vertical grid points

    # Load B-field components and replace NaNs
    bz = np.nan_to_num(nc.variables['bz'][:]).reshape(len(R), len(z))
    bt = np.nan_to_num(nc.variables['bt'][:]).reshape(len(R), len(z))
    br = np.nan_to_num(nc.variables['br'][:]).reshape(len(R), len(z))

    # Load densities and reshape
    ne = np.nan_to_num(nc.variables['ne'][:]).reshape(len(R), len(z))
    ni = np.nan_to_num(nc.variables['ni'][:]).reshape(len(R), len(z))

    # Load temperatures and reshape
    te = np.nan_to_num(nc.variables['te'][:]).reshape(len(R), len(z))
    ti = np.nan_to_num(nc.variables['ti'][:]).reshape(len(R), len(z))

    # Load velocities and reshape
    vt = np.nan_to_num(nc.variables['vt'][:]).reshape(len(R), len(z))
    vr = np.nan_to_num(nc.variables['vr'][:]).reshape(len(R), len(z))
    vz = np.nan_to_num(nc.variables['vz'][:]).reshape(len(R), len(z))

# Load centroid and norm vector data
centroid = np.loadtxt('centroid.csv', delimiter=',')
r_centroid = np.sqrt(centroid[:, 0]**2 + centroid[:, 1]**2)
norm_vec = np.loadtxt('norm_vec.csv', delimiter=',')

# Define interpolators for profiles
interp_ne = RegularGridInterpolator((R, z), ne)
interp_ni = RegularGridInterpolator((R, z), ni)
interp_te = RegularGridInterpolator((R, z), te)
interp_ti = RegularGridInterpolator((R, z), ti)
interp_br = RegularGridInterpolator((R, z), br)
interp_bt = RegularGridInterpolator((R, z), bt)
interp_bz = RegularGridInterpolator((R, z), bz)
interp_vr = RegularGridInterpolator((R, z), vr)
interp_vt = RegularGridInterpolator((R, z), vt)
interp_vz = RegularGridInterpolator((R, z), vz)

# Interpolate profiles onto the surface
points = np.vstack((r_centroid, centroid[:, 2])).T
ne_surf = interp_ne(points)
ni_surf = interp_ni(points)
te_surf = interp_te(points)
ti_surf = interp_ti(points)
br_surf = interp_br(points)
bt_surf = interp_bt(points)
bz_surf = interp_bz(points)
vr_surf = interp_vr(points)
vt_surf = interp_vt(points)
vz_surf = interp_vz(points)

# Calculate B-field components in Cartesian coordinates
phi_centroid = np.arctan2(centroid[:, 1], centroid[:, 0])
bx = br_surf * np.cos(phi_centroid) - bt_surf * np.sin(phi_centroid)
by = br_surf * np.sin(phi_centroid) + bt_surf * np.cos(phi_centroid)
bz = bz_surf
b_mag = np.sqrt(bx**2 + by**2 + bz**2)

# Normalize B-field components
ubx = bx / b_mag
uby = by / b_mag
ubz = bz / b_mag

# Calculate velocity components in Cartesian coordinates
vx = vr_surf * np.cos(phi_centroid) - vt_surf * np.sin(phi_centroid)
vy = vr_surf * np.sin(phi_centroid) + vt_surf * np.cos(phi_centroid)
vz = vz_surf
v_mag = np.sqrt(vx**2 + vy**2 + vz**2)

# Normalize velocity components
uvx = vx / v_mag
uvy = vy / v_mag
uvz = vz / v_mag

# Calculate flux and save it
o8plus_flux_surf = ni_surf * v_mag
np.savetxt('o8plus_flux_surf.csv', o8plus_flux_surf, delimiter=',')

# Calculate theta
norm_vec_mag = np.linalg.norm(norm_vec, axis=1)
unorm_vec = norm_vec / norm_vec_mag[:, None]
theta = np.arccos(unorm_vec[:, 0] * ubx + unorm_vec[:, 1] * uby + unorm_vec[:, 2] * ubz)
theta = np.where(theta > np.pi / 2, np.abs(theta - np.pi), theta)

# Save selected profiles to CSV
np.savetxt('ne_surf.csv', ne_surf, delimiter=',')
np.savetxt('te_surf.csv', te_surf, delimiter=',')
np.savetxt('ti_surf.csv', ti_surf, delimiter=',')
np.savetxt('theta.csv', theta, delimiter=',')



# Histogram of theta
plt.figure()
plt.hist(theta, bins=50)
plt.xlabel("Theta (radians)")
plt.ylabel("Frequency")
plt.title("Histogram of Theta")
plt.show()


# Sheath Calculations
sheathType = 1
potential_surf = np.zeros(len(centroid))
if sheathType == 1:  # Thermal sheath
    me = 1 / 2000
    background_amu = 2
    sheath_factor = np.abs(0.5 * np.log((2 * np.pi * me / background_amu) * (1 + ti_surf / te_surf)))
    potential_surf = sheath_factor * ti_surf

# Write surface data
surface_data = {
    "potential_surf": potential_surf,
    "ne_surf": ne_surf,
    "te_surf": te_surf,
    "ti_surf": ti_surf,
    "vp_surf": vz_surf,
    "b_mag": b_mag,
    "br_surf": br_surf,
    "bt_surf": bt_surf,
    "bz_surf": bz_surf,
    "theta": r_centroid,
    "ni_surf": ni_surf
}
surface_variables = pd.DataFrame(surface_data)
surface_variables.to_csv('Targets.txt', sep='\t', index=False)

print("Script execution complete.")


# Generate GITR geometry configuration file
with open('gitrGeometryPointPlane3d.cfg', 'w', encoding='utf-8') as f:
    f.write("geom = {\n   x1 = [")
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 0]))
    f.write("]\n   y1 = [")
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 1]))
    f.write("]\n   z1 = [")
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 2]))
    # Repeat similar blocks for x2, y2, z2, x3, y3, z3, etc.
    f.write("]\n   inDir = [")
    f.write(', '.join(f'{int(x)}' for x in inDir))
    f.write("]\n   potential = [")
    f.write(', '.join(f'{0.0:.5e}' for _ in range(len(planes))))
    f.write("]\n}")


import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from netCDF4 import Dataset

import numpy as np
from netCDF4 import Dataset

print(">>>> Initializing Particles")

# Assuming `surfs`, `vp_surf`, `area`, `X`, `Y`, `Z`, `abcd`, `plane_norm`, and `inDir` are defined

surf_inds = np.where(surfs)[0]
nP = 1000
nR = 100

# Erosion data
Y0 = 0.1
eroded_flux = Y0 * vz_surf  # Update with realistic values
erosion = eroded_flux * area
erosion = erosion[surf_inds]
erosion_rate = np.sum(erosion)

# Cumulative distribution function for erosion
erosion_inds = np.where(erosion > 0)[0]
erosion_sub = erosion[erosion_inds]
erosion_sub_cdf = np.cumsum(erosion_sub)
erosion_rate = erosion_sub_cdf[-1]
erosion_sub_cdf /= erosion_sub_cdf[-1]

# Random selection of elements
rand1 = np.random.rand(nP)
element = np.interp(rand1, np.hstack(([0], erosion_sub_cdf)), np.arange(len(erosion_sub_cdf)))
element_ceil = np.ceil(element).astype(int) - 1

# Particle initialization
xP = np.zeros(nP)
yP = np.zeros(nP)
zP = np.zeros(nP)
vxP = np.zeros(nP)
vyP = np.zeros(nP)
vzP = np.zeros(nP)

offset = 1e-5

def sample_triangle(x, y, z, n_samples):
    # Randomly sample points inside a triangle
    a = np.random.rand(n_samples, 1)
    b = np.random.rand(n_samples, 1)
    mask = a + b > 1
    a[mask] = 1 - a[mask]
    b[mask] = 1 - b[mask]
    return (1 - a - b) * x[0] + a * x[1] + b * x[2], \
           (1 - a - b) * y[0] + a * y[1] + b * y[2], \
           (1 - a - b) * z[0] + a * z[1] + b * z[2]

# Calculate particle positions and velocities
for j in range(nP):
    i = erosion_inds[element_ceil[j]]
    normal = -abcd[i, :3] / plane_norm[i]
    normal *= inDir[i]

    # Apply offset along normal
    x_tri = X[i, :] + offset * normal[0]
    y_tri = Y[i, :] + offset * normal[1]
    z_tri = Z[i, :] + offset * normal[2]

    # Sample point inside triangle
    xP[j], yP[j], zP[j] = sample_triangle(x_tri, y_tri, z_tri, 1)
    vxP[j] = 5000 * normal[0]
    vyP[j] = 5000 * normal[1]
    vzP[j] = 5000 * normal[2]

# Visualization (optional)
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X.flatten(), Y.flatten(), Z.flatten(), color='b', alpha=0.3)
ax.scatter(xP, yP, zP, color='r')
plt.show()

# Saving to NetCDF file
with Dataset('particle_source_west.nc', 'w', format='NETCDF4') as ncid:
    ncid.createDimension('nP', nP)

    xVar = ncid.createVariable('x', 'f8', ('nP',))
    yVar = ncid.createVariable('y', 'f8', ('nP',))
    zVar = ncid.createVariable('z', 'f8', ('nP',))
    vxVar = ncid.createVariable('vx', 'f8', ('nP',))
    vyVar = ncid.createVariable('vy', 'f8', ('nP',))
    vzVar = ncid.createVariable('vz', 'f8', ('nP',))

    xVar[:] = xP
    yVar[:] = yP
    zVar[:] = zP
    vxVar[:] = vxP
    vyVar[:] = vyP
    vzVar[:] = vzP

print("Particle source saved to 'particle_source_west.nc'")
