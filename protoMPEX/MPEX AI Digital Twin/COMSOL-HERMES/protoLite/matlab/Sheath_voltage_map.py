# plot_sheath_voltage_mpex.py
# Computes and plots V_RF and V_DC from Efield.csv and window_parameters.csv
# Consistent with Sheath_impedance_MPEXAI.py

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, atan2, tanh, log
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

# ----------------------------
# User settings
# ----------------------------
efield_file = "Efield.csv"
window_file = "window_parameters.csv"

Te = 8.0          # eV
d = 0.002         # m, RF layer thickness
smooth_sigma = 2.0

# Set to True only if row order is identical between files
use_row_by_row_pairing = True

# ----------------------------
# Read files
# ----------------------------
def read_complex_csv(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            row = []
            for v in line.split(","):
                v = v.strip()
                if "i" in v:
                    row.append(complex(v.replace("i", "j")))
                else:
                    row.append(float(v))
            data.append(row)
    return np.array(data, dtype=complex)

def read_real_csv(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            data.append([float(v.strip()) for v in line.split(",")])
    return np.array(data, dtype=float)

LayerData = read_complex_csv(efield_file)
EdgeData = read_real_csv(window_file)

# Efield columns
x = LayerData[:, 0].real
y = LayerData[:, 1].real
z = LayerData[:, 2].real
Ex = LayerData[:, 3]
Ey = LayerData[:, 4]

if use_row_by_row_pairing:
    N = min(len(LayerData), len(EdgeData))
    x = x[:N]
    y = y[:N]
    z = z[:N]
    Ex = Ex[:N]
    Ey = Ey[:N]

    Density = EdgeData[:N, 3]
    Br = EdgeData[:N, 4]
    B0 = EdgeData[:N, 5]
else:
    # Interpolate window data onto Efield coordinates
    xw = EdgeData[:, 0]
    yw = EdgeData[:, 1]
    zw = EdgeData[:, 2]

    Density = griddata((xw, yw, zw), EdgeData[:, 3], (x, y, z), method="nearest")
    Br = griddata((xw, yw, zw), EdgeData[:, 4], (x, y, z), method="nearest")
    B0 = griddata((xw, yw, zw), EdgeData[:, 5], (x, y, z), method="nearest")

# ----------------------------
# Geometry and V_RF
# ----------------------------
r = np.sqrt(x**2 + y**2)
phi = np.arctan2(y, x)
phi_deg = phi * 180.0 / np.pi

E_norm = (Ex * x + Ey * y) / r
V_RF = np.abs(E_norm) * d

# ----------------------------
# V_DC model from Sheath_impedance_MPEXAI.py
# ----------------------------
def compute_phi0avg(Vlayer, Density, Br, B0, Te):
    mu = 24.17
    Z = 1
    A = 2
    w_rf = 13.56e6 * 2.0 * pi

    omegapi = 1.32e3 * Z * np.sqrt(Density * 1e-6 / A)
    omega_hat = w_rf / omegapi

    xi = Vlayer / Te

    j = 0.0
    upar0 = 1.1

    a1 = 3.70285
    a2 = 3.81991
    b1 = 1.13352
    b2 = 1.24171
    a3 = 2.0 * b2 / pi

    c0 = 0.966463
    c1 = 0.141639

    gg = c0 + c1 * np.tanh(omega_hat)
    xi1 = gg * xi

    phi0avg = (
        (np.log(mu) + xi1 * a1 + xi1**2 * a2 + xi1**3 * a3)
        / (1.0 + xi1 * b1 + xi1**2 * b2)
        - np.log(1.0 - j / upar0)
        + np.log(mu / 24.17)
    )

    return np.real(phi0avg)

phi0avg = compute_phi0avg(V_RF, Density, Br, B0, Te)
V_DC = phi0avg * Te
ratio = V_RF / V_DC

# ----------------------------
# Save data
# ----------------------------
out = np.column_stack((phi_deg, z, V_RF, V_DC, ratio, Density, Br, B0))
np.savetxt(
    "sheath_voltage_phi_z_python.csv",
    out,
    delimiter=",",
    header="phi_deg,z,V_RF,V_DC,V_RF_over_V_DC,Density,Br,B0",
    comments="",
)

print("Saved sheath_voltage_phi_z_python.csv")
print(f"V_RF range: {np.nanmin(V_RF):.3g} to {np.nanmax(V_RF):.3g} V")
print(f"V_DC range: {np.nanmin(V_DC):.3g} to {np.nanmax(V_DC):.3g} V")
print(f"V_RF/V_DC range: {np.nanmin(ratio):.3g} to {np.nanmax(ratio):.3g}")

# ----------------------------
# Plot helper
# ----------------------------
def plot_map(phi_deg, z, val, title, cbar_label, fname, vmin=None, vmax=None):
    good = np.isfinite(phi_deg) & np.isfinite(z) & np.isfinite(val)

    p = phi_deg[good]
    zz = z[good]
    vv = val[good]

    p_grid = np.linspace(-180, 180, 500)
    z_grid = np.linspace(np.min(zz), np.max(zz), 350)
    P, Z = np.meshgrid(p_grid, z_grid)

    V = griddata((p, zz), vv, (P, Z), method="linear")
    Vn = griddata((p, zz), vv, (P, Z), method="nearest")
    V[np.isnan(V)] = Vn[np.isnan(V)]

    if smooth_sigma > 0:
        V = gaussian_filter(V, smooth_sigma)

    plt.figure(figsize=(10, 7))
    im = plt.imshow(
        V,
        origin="lower",
        extent=[p_grid.min(), p_grid.max(), z_grid.min(), z_grid.max()],
        aspect="auto",
        cmap="turbo",
        vmin=vmin,
        vmax=vmax,
    )

    plt.colorbar(im, label=cbar_label)
    plt.xlabel(r"Azimuthal coordinate $\phi$ [deg]", fontsize=15, fontweight="bold")
    plt.ylabel("Axial coordinate z [m]", fontsize=15, fontweight="bold")
    plt.title(title, fontsize=17, fontweight="bold")
    plt.xlim(-180, 180)
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.show()

# ----------------------------
# Make plots
# ----------------------------
plot_map(
    phi_deg,
    z,
    V_RF,
    r"RF Sheath Voltage $V_{RF}(\phi,z)$",
    r"$V_{RF}$ [V]",
    "V_RF_phi_z_python.png",
)

plot_map(
    phi_deg,
    z,
    V_DC,
    r"DC Sheath Voltage $V_{DC}(\phi,z)$",
    r"$V_{DC}$ [V]",
    "V_DC_phi_z_python.png",
)

plot_map(
    phi_deg,
    z,
    ratio,
    r"$V_{RF}/V_{DC}(\phi,z)$",
    r"$V_{RF}/V_{DC}$",
    "V_RF_over_V_DC_phi_z_python.png",
)

# ----------------------------
# Histograms
# ----------------------------
plt.figure(figsize=(8, 5))
plt.hist(V_DC[np.isfinite(V_DC)], bins=60)
plt.xlabel(r"$V_{DC}$ [V]", fontsize=14, fontweight="bold")
plt.ylabel("Counts", fontsize=14, fontweight="bold")
plt.title(r"Distribution of $V_{DC}$", fontsize=16, fontweight="bold")
plt.tight_layout()
plt.savefig("V_DC_hist_python.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(V_RF[np.isfinite(V_RF)], bins=60)
plt.xlabel(r"$V_{RF}$ [V]", fontsize=14, fontweight="bold")
plt.ylabel("Counts", fontsize=14, fontweight="bold")
plt.title(r"Distribution of $V_{RF}$", fontsize=16, fontweight="bold")
plt.tight_layout()
plt.savefig("V_RF_hist_python.png", dpi=300)
plt.show()