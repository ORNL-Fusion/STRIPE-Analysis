import os
import subprocess
import impacts
import shutil
import numpy as np
target = np.loadtxt('Targets.txt', dtype='float',skiprows=1,delimiter=',')
print(target,target.shape)

for i in range(0, 14583):
    for j in range(8,9):
        surface_potential=target[i,0]
        ne = target[i,1]
        te = target[i,2]
        ti = target[i,3]
        vp = target[i,4]
        btot = target[i,5]
        br = target[i,6]
        bt = target[i,7]
        bz = target[i,8]
        theta = target[i,9]
        #impacts.d3d_case_C(charge = j, Tee=te, Tii = ti, n=ne,btot = btot,br = br, bt = bt, bz = bz, sheath_factor = 1.0)
        impacts.west_thermal_O(charge = 8,Tee=te, Tii=ti, n=ne,btot=btot,theta=theta, sheath_factor=1.0)
        #impacts.west_rf_O(surface_potential = surface_potential, charge = 1, Tee=te, Tii = ti, n=ne,btot = btot,theta=theta, sheath_factor = 1.0)
        #impacts.iter_case_D(surface_potential = surface_potential, charge = 1, Tee=te, Tii = ti, n=ne,btot = btot,theta=theta, sheath_factor = 1.0)
        #subprocess.run("/global/cfs/cdirs/m77/atul/gitr/build/./GITR", shell=True, check=True)        
        subprocess.run("~/myRepos/gitr_new/build/GITR", shell=True, check=True)
        shutil.copyfile("output/surface.nc","surface_C"+str(j)+"_loc_"+str(i) +".nc")
        shutil.copyfile("output/positions.nc","positions_C"+str(j)+"_loc_"+str(i)+".nc")
        os.remove("output/surface.nc")
        os.remove("output/positions.nc")
        os.remove("output/positions.m")
        #os.remove("output/particleSource.nc")
