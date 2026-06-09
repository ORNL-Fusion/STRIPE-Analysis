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



for filename in files:
    T = Torosurf.loadtxt(os.path.join("../sealed", filename))
    T.phi += 180
    T.rz = np.flip(T.rz, 1)
    T.savetxt(filename, fmt='%11.3f', skipv=True)
