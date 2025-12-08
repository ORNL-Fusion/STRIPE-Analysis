import pickle
import numpy as np
import matplotlib.pyplot as plt

# — Load the pickle —
with open("profs_200882_1800_py3.pkl", "rb") as f:
    data = pickle.load(f)

# — Raw data (nedatpsi) —
ned     = data["nedatpsi"]
psi_raw = np.array(ned["x"])
ne_raw  = np.array(ned["y"])
ne_err  = np.array(ned.get("err", np.full_like(ne_raw, 0.05)))

# — Fitted profile (netanhpsi) —
tanh    = data["netanhpsi"]
psi_f   = np.array(tanh["x"])
ne_f    = np.array(tanh["y"])

# — Make two stacked subplots — 
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 10), sharex=False)

# 1) Raw data on its own ψ-grid
ax1.errorbar(
    psi_raw, ne_raw, yerr=ne_err,
    fmt="o", capsize=3, label="Raw $n_e$"
)
ax1.set_xlabel("Normalized poloidal flux $\\psi_N$")
ax1.set_ylabel("Electron density $n_e$ [$10^{20}\\,$m$^{-3}$]")
ax1.set_title("Raw profile (nedatpsi)")
ax1.grid(True)
ax1.legend()

# 2) Fitted profile on its own ψ-grid
ax2.plot(
    psi_f, ne_f,
    "r--", lw=2, label="Fit $n_e$"
)
ax2.set_xlabel("Normalized poloidal flux $\\psi_N$")
ax2.set_ylabel("Electron density $n_e$ [$10^{20}\\,$m$^{-3}$]")
ax2.set_title("Tanh‐fit profile (netanhpsi)")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
