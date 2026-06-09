import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

import sys

if len(sys.argv) != 2:
    print(f"Usage:\n  {sys.argv[0]} FILE.nc")

filename = sys.argv[1]

with nc.Dataset(filename, mode='r') as f:
    r = f.variables['r'][:]
    z = f.variables['z'][:]
    angle = f.variables['angle'][:]
    Ne = f.variables['electron_density'][:] # m^-3
    Te = f.variables['electron_temperature'][:]  # eV

# Make Cartesian coordinates
X = r * np.cos(angle)
Y = r * np.sin(angle)

yind = 20  # index along the axis (Y)

fig, ax = plt.subplots()
cm = ax.contourf(X[:,yind,:], Y[:,yind,:], Ne[:,yind,:], 50)
plt.colorbar(cm, ax=ax)
ax.set_aspect('equal')
ax.set_title(f"Plasma density at Z={z[0,yind,0]:.2f}m")
ax.set_xlabel("[m]")
ax.set_ylabel("[m]")
fig.savefig("density.png")
fig.savefig("density.pdf")
plt.show()
plt.close(fig)

fig, ax = plt.subplots()
cm = ax.contourf(X[:,yind,:], Y[:,yind,:], Te[:,yind,:], 50)
plt.colorbar(cm, ax=ax)
ax.set_aspect('equal')
ax.set_title(f"Electron temperature [eV] at Z={z[0,yind,0]:.2f}m")
ax.set_xlabel("[m]")
ax.set_ylabel("[m]")
fig.savefig("temperature.png")
fig.savefig("temperature.pdf")
plt.show()
plt.close(fig)
