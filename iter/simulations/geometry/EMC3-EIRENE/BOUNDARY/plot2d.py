import numpy as np
import matplotlib.pyplot as plt
import os.path

from moose.geometry import Torosurf



files = [
    'CENTRALLOWERTHIRD_MIRR_ET.tgt',
    'CENTRALMEDTHIRD_MIRR.tgt',
    'CENTRALUPPERTHIRD_MIRR_ET.tgt',
    'FULLLEFT_ET.tgt',
    'FULLLEFT_MIRR_ET.tgt',
    'LOWERHALFLEFT_ET.tgt',
    'LOWERHALFLEFT_MIRR_ET.tgt',
    'UPPERHALFLEFT_ET.tgt',
    'UPPERHALFLEFT_MIRR_ET.tgt',
    'CENTRALMIDSEAL.txt',
    'CENTRALUPSEAL.txt'
    ]



T = {filename: Torosurf.loadtxt(filename) for filename in files}



def plot2d(T, iphi):
    #for i in range(T.nu):
        #phi = T.phi[i]
        #plt.plot(T.r[:,i], T.z[:,i], label=T.label+f", phi = {phi}")
    plt.plot(T.r[:,iphi], T.z[:,iphi])
    print(T.phi[iphi])



plot2d(T["LOWERHALFLEFT_ET.tgt"],-1)
plot2d(T['CENTRALLOWERTHIRD_MIRR_ET.tgt'], 0)
#plt.legend()
plt.show()
