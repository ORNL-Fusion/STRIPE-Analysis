import numpy as np
import os
from scipy.interpolate import interpn

# Read in the data from MATLAB files
planes = np.genfromtxt('west_geom.mat', delimiter=',')
centroid = np.loadtxt('centroid.csv')
norm_vec = np.loadtxt('norm_vec.csv')

# Calculate surface properties
inDir = np.ones(len(planes))
norm_vec_mag = np.sqrt(norm_vec[:, 0]**2 + norm_vec[:, 1]**2 + norm_vec[:, 2]**2)
unorm_vec = norm_vec / norm_vec_mag[:, None]

for i in range(len(planes)):
    normal = -planes[i, 0:3] / np.linalg.norm(planes[i, 0:3])
    dot_product = np.dot(normal, centroid[i])
    if dot_product > 0:
        inDir[i] = -1
    normal *= inDir[i]
    
    if abs(normal[2]) > 0.005:
        planes[i, -1] = 0
        planes[i, -2] = 0

np.savetxt('centroid.csv', centroid, delimiter=',')
np.savetxt('norm_vec.csv', norm_vec, delimiter=',')

exit()

# Generate GITR geometry configuration file
with open('gitrGeometryPointPlane3d.cfg', 'w') as f:
    f.write('geom = \n{ \n   x1 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 0]))
    f.write('] \n   y1 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 1]))
    f.write('] \n   z1 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 2]))
    f.write('] \n   x2 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 3]))
    f.write('] \n   y2 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 4]))
    f.write('] \n   z2 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 5]))
    f.write('] \n   x3 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 6]))
    f.write('] \n   y3 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 7]))
    f.write('] \n   z3 = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 8]))
    f.write('] \n   a = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 0]))
    f.write('] \n   b = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 1]))
    f.write('] \n   c = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 2]))
    f.write('] \n   d = [')
    f.write(', '.join(f'{x:.5e}' for x in planes[:, 3]))
    f.write('] \n   plane_norm = [')
    f.write(', '.join(f'{x:.5e}' for x in np.linalg.norm(planes[:, 0:3], axis=1)))
    f.write('] \n   BCxBA = [')
    f.write(', '.join(f'{x:.5e}' for x in np.cross(planes[:, 3:6] - planes[:, 0:3], planes[:, 6:9] - planes[:, 0:3], axis=1)))
    f.write('] \n   CAxCB = [')
    f.write(', '.join(f'{x:.5e}' for x in np.cross(planes[:, 6:9] - planes[:, 3:6], planes[:, 0:3] - planes[:, 3:6], axis=1)))
    f.write('] \n   area = [')
    f.write(', '.join(f'{x:.5e}' for x in 0.5 * np.linalg.norm(np.cross(planes[:, 3:6] - planes[:, 0:3], planes[:, 6:9] - planes[:, 0:3], axis=1), axis=1)))
    f.write('] \n   Z = [')
    f.write(', '.join(f'{x:.1f}' for x in planes[:, -2]))
    f.write('] \n   surface = [')
    f.write(', '.join(f'{int(x)}' for x in planes[:, -1]))
    f.write('] \n   inDir = [')
    f.write(', '.join(f'{int(x)}' for x in inDir))
    f.write('] \n   potential = [')
    f.write(', '.join(f'{0.0:.5e}' for _ in range(len(planes)))) 
    f.write('] \n}')
    f.write('\nperiodic = 0;\ntheta0 = 0.0;\ntheta1 = 0.0\nperiodic_bc_x0 = 0.0;\nperiodic_bc_x1 = 0.0;\nperiodic_bc_x = 0;}')

# Generate particle initialization
nP = 1000
nR = 100
Y0 = 0.1
eroded_flux = Y0 * np.loadtxt('o8plus_flux_surf.csv')
erosion = eroded_flux * 0.5 * np.linalg.norm(np.cross(planes[:, 3:6] - planes[:, 0:3], planes[:, 6:9] - planes[:, 0:3], axis=1), axis=1)
erosion_rate = np.sum(erosion)

erosion_inds = np.nonzero(erosion)[0]
erosion_sub = erosion[erosion_inds]
erosion_sub_cdf = np.cumsum(erosion_sub) / np.sum(erosion_sub)

rand1 = np.random.rand(nP)
element = np.interp(rand1, [0, 1], [0, len(erosion_sub_cdf) - 1])
element_ceil = np.ceil(element).astype(int)

xP = np.zeros(nP)
yP = np.zeros(nP)
zP = np.zeros(nP)
vxP = np.zeros(nP)
vyP = np.zeros(nP)
vzP = np.zeros(nP)

offset = 1e-5
for j in range(nP):
    i = erosion_inds[element_ceil[j] - 1]
    normal = -planes[i, 0:3] / np.linalg.norm(planes[i, 0:3])
    normal *= inDir[i]

    x_tri = planes[i, 0:3] + offset * normal
    y_tri = planes[i, 3:6] + offset * normal
    z_tri = planes[i, 6:9] + offset * normal

    samples = sample_triangle(x_tri, y_tri, z_tri, 1)
    xP[j] = samples[0, 0]
    yP[j] = samples[0, 1]
    zP[j] = samples[0, 2]
    vxP[j] = 5000 * normal[0]
    vyP[j] = 5000 * normal[1] 
    vzP[j] = 5000 * normal[2]

np.savez('particle_source_west.nc', x=xP, y=yP, z=zP, vx=vxP, vy=vyP, vz=vzP)


def sample_triangle(x, y, z, nP):
    x_transform = x - x[0]
    y_transform = y - y[0]
    z_transform = z - z[0]

    v1 = [x_transform[1], y_transform[1], z_transform[1]]
    v2 = [x_transform[2], y_transform[2], z_transform[2]]
    v12 = np.array(v2) - np.array(v1)
    normal = np.cross(v1, v2)

    a1 = np.random.rand(nP)
    a2 = np.random.rand(nP)

    samples = a1[:, None] * np.array(v1) + a2[:, None] * np.array(v2)
    samples2x = samples[:, 0] - v2[0]
    samples2y = samples[:, 1] - v2[1]
    samples2z = samples[:, 2] - v2[2]
    samples12x = samples[:, 0] - v1[0]
    samples12y = samples[:, 1] - v1[1]
    samples12z = samples[:, 2] - v1[2]
    v1Cross = [(v1[1] * samples[:, 2] - v1[2] * samples[:, 1]), 
               (v1[2] * samples[:, 0] - v1[0] * samples[:, 2]),
               (v1[0] * samples[:, 1] - v1[1] * samples[:, 0])]
    v2 = -np.array(v2)
    v2Cross = [(v2[1] * samples2z - v2[2] * samples2y),
               (v2[2] * samples2x - v2[0] * samples2z),
               (v2[0] * samples2y - v2[1] * samples2x)]
    v12Cross = [(v12[1] * samples12z - v12[2] * samples12y),
                (v12[2] * samples12x - v12[0] * samples12z),
                (v12[0] * samples12y - v12[1] * samples12x)]

    v1CD = normal[0] * v1Cross[0] + normal[1] * v1Cross[1] + normal[2] * v1Cross[2]
    v2CD = normal[0] * v2Cross[0] + normal[1] * v2Cross[1] + normal[2] * v2Cross[2]
    v12CD = normal[0] * v12Cross[0] + normal[1] * v12Cross[1] + normal[2] * v12Cross[2]

    inside = np.abs(np.sign(v1CD) + np.sign(v2CD) + np.sign(v12CD))
    insideInd = np.where(inside == 3)[0]
    notInsideInd = np.where(inside != 3)[0]

    v2 = -np.array(v2)
    dAlongV1 = v1[0] * samples[notInsideInd, 0] + v1[1] * samples[notInsideInd, 1] + v1[2] * samples[notInsideInd, 2]
    dAlongV2 = v2[0] * samples[notInsideInd, 0] + v2[1] * samples[notInsideInd, 1] + v2[2] * samples[notInsideInd, 2]

    dV1 = np.linalg.norm(v1)
    dV2 = np.linalg.norm(v2)
    halfdV1 = 0.5 * dV1
    halfdV2 = 0.5 * dV2

    samples[notInsideInd] = [-(samples[notInsideInd, 0] - 0.5 * v1[0]) + 0.5 * v1[0],
                             -(samples[notInsideInd, 1] - 0.5 * v1[1]) + 0.5 * v1[1],
                             -(samples[notInsideInd, 2] - 0.5 * v1[2]) + 0.5 * v1[2]]
    samples[notInsideInd] = [(samples[notInsideInd, 0] + v2[0]),
                             (samples[notInsideInd, 1] + v2[1]),
                             (samples[notInsideInd, 2] + v2[2])]

    samples[:, 0] += x[0]
    samples[:, 1] += y[0]
    samples[:, 2] += z[0]

    return samples