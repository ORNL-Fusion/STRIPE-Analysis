import matplotlib.pyplot as plt
import numpy as np
import adas_parser

data = adas_parser.load('pec93#o_pju#o1', return_xr=True)
for i, ne in enumerate([1e17, 1e20]):
    da = data.sel(line='4419.3 A').sel(ne=ne, method='nearest')
    da.sel(TYPE='EXCIT').plot(label='EXCIT: ne={}'.format(da['ne'].item()), color='C{}'.format(i))
    da.sel(TYPE='RECOM').plot(label='RECOM: ne={}'.format(da['ne'].item()), color='C{}'.format(i), ls='--')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.ylim(1e-22, 1e-14)