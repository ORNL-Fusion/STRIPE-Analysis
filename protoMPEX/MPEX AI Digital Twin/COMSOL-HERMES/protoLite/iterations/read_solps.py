from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


rootgrp = Dataset("protoMPEX_profiles.nc", "r", format="NETCDF4")
#print(rootgrp.variables)

fig, ax = plt.subplots()
pcm=ax.pcolormesh(rootgrp.variables['r'],rootgrp.variables['z'],rootgrp.variables['ne']) #,norm=LogNorm()
cbar = fig.colorbar(pcm, ax=ax)
plt.show()

rkron=np.kron(np.ones(rootgrp.variables['z'].shape),rootgrp.variables['r']);
zkron=np.kron(rootgrp.variables['z'],np.ones(rootgrp.variables['r'].shape))
np.savetxt('ne_from_nc.txt',np.stack((rkron,zkron,np.array(rootgrp.variables['ne']).flatten(),np.array(rootgrp.variables['br']).flatten(),np.array(rootgrp.variables['bt']).flatten(),np.array(rootgrp.variables['bz']).flatten(),np.array(rootgrp.variables['te']).flatten()),axis=-1))

rootgrp.close()

#rootgrp = Dataset("bfield_protoMPEX.nc", "r", format="NETCDF4")
#print(rootgrp.variables)
#rootgrp.close()
exit()
