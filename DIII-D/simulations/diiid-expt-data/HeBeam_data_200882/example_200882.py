from uncertainties import unumpy as unp
import numpy as np
import matplotlib.pyplot as plt

data = np.load("He-Beam_200882.npy", allow_pickle=True).item()
quickfit = np.load("quickfit_200882.npy", allow_pickle=True).item()


plt.pcolormesh(data["R"], data["time"], data["T_e"].T)
plt.colorbar(label="$T_e$ (eV)")
plt.ylim(1000,4000)
plt.ylabel("time (ms)")
plt.xlabel("R (m)")
plt.show()

plt.pcolormesh(data["R"], data["time"], data["n_e"].T)
plt.colorbar(label="$n_e$ (m$^{-3}$)")
plt.ylim(1000,4000)
plt.ylabel("time (ms)")
plt.xlabel("R (m)")
plt.show()


time_slice = 2600 
loc_quickfit = np.argmin(np.abs(quickfit["time"]-time_slice))
loc_hebeam = np.argmin(np.abs(data["time"]-time_slice))

fig, ax = plt.subplots(2, 1, figsize=(6,8), sharex=True)
ax[0].plot(quickfit["psi_n"][loc_quickfit], unp.nominal_values(quickfit["T_e"][loc_quickfit]))
ax[0].fill_between(quickfit["psi_n"][loc_quickfit], unp.nominal_values(quickfit["T_e"][loc_quickfit])- unp.std_devs(quickfit["T_e"][loc_quickfit]), unp.nominal_values(quickfit["T_e"][loc_quickfit])+unp.std_devs(quickfit["T_e"][loc_quickfit]), alpha=0.3)
ax[0].set_ylabel("$T_e$ (eV)")

ax[1].plot(quickfit["psi_n"][loc_quickfit], unp.nominal_values(quickfit["n_e"][loc_quickfit]))
ax[1].fill_between(quickfit["psi_n"][loc_quickfit], unp.nominal_values(quickfit["n_e"][loc_quickfit])- unp.std_devs(quickfit["n_e"][loc_quickfit]), unp.nominal_values(quickfit["n_e"][loc_quickfit])+unp.std_devs(quickfit["n_e"][loc_quickfit]), alpha=0.3)
ax[1].set_xlabel("$\psi_n$")
ax[1].set_ylabel("$n_e$ (eV)")
plt.title(f"200882@{time_slice} ms")
plt.show()


fig, ax = plt.subplots(2, 1, figsize=(6,8), sharex=True)
ax[0].plot(data["psi_n"][:,loc_hebeam], data["T_e"][:,loc_hebeam])
ax[0].set_ylabel("$T_e$ (eV)")

ax[1].plot(data["psi_n"][:,loc_hebeam], data["n_e"][:,loc_hebeam])
ax[1].set_xlabel("$\psi_n$")
ax[1].set_ylabel("$n_e$ (eV)")
plt.title(f"200882@{time_slice} ms")
plt.show()