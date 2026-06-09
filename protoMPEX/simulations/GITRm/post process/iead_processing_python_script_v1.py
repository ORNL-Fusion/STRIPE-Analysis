#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LogNorm

# ============================
# CONFIG
# ============================
tilt = 85
file = f"../tilted_targets/test/{tilt}_degrees/MPEX_runs_v1/gitrm-surface.nc"

SKIP_FIRST_NSURFACES = 6
DROP_LAST_NSURFACES = 3

COARSEN_ENERGY_FACTOR = 5
COARSEN_ANGLE_FACTOR = 1

scale = 1.0

BIN_COMBINED = 1
BIN_TA = [3, 4, 5, 6]

E_MIN = 0.0
E_MAX = 500.0

A_MIN = 0.0
A_MAX = 90.0

SMOOTH_SIGMA = 0.8

# Log-scale color limits
CLIM_INC_LOG = [1e-5, 1e1]
CLIM_REFL_LOG = [1e-8, 1e-2]

X_LIMIT = [0, 100]
Y_LIMIT = [0, 90]

# ============================
# HELPERS
# ============================
def block_mean_1d(x, f):
    x = np.asarray(x)
    n2 = (len(x) // f) * f
    x = x[:n2]
    return x.reshape(-1, f).mean(axis=1)


def block_mean_2d(Z, fa, fe):
    nA2 = (Z.shape[0] // fa) * fa
    nE2 = (Z.shape[1] // fe) * fe
    Z = Z[:nA2, :nE2]
    Z = Z.reshape(nA2 // fa, fa, nE2 // fe, fe)
    return np.nanmean(np.nanmean(Z, axis=3), axis=1)


def prep_AES_to_AE(AES):
    # AES expected as [angle, energy, surface]
    nA, nE, nS = AES.shape

    i0 = SKIP_FIRST_NSURFACES
    i1 = nS - DROP_LAST_NSURFACES

    if i0 >= i1:
        raise ValueError("Invalid surface selection.")

    AE = np.nansum(AES[:, :, i0:i1], axis=2)

    angle_deg = np.linspace(0, 90, nA)
    energy_eV = np.linspace(0, 500, nE)

    fa = max(1, int(round(COARSEN_ANGLE_FACTOR)))
    fe = max(1, int(round(COARSEN_ENERGY_FACTOR)))

    nA2 = (AE.shape[0] // fa) * fa
    nE2 = (AE.shape[1] // fe) * fe

    AE = AE[:nA2, :nE2]
    angle_deg = angle_deg[:nA2]
    energy_eV = energy_eV[:nE2]

    AE = block_mean_2d(AE, fa, fe)
    angle_deg = block_mean_1d(angle_deg, fa)
    energy_eV = block_mean_1d(energy_eV, fe)

    AE = scale * AE
    AE[AE <= 0] = np.nan

    return AE, angle_deg, energy_eV


def smooth_for_plot(Z):
    Zfill = np.array(Z, dtype=float)
    Zfill[~np.isfinite(Zfill)] = 0.0

    Zsmooth = gaussian_filter(Zfill, sigma=SMOOTH_SIGMA)

    Zsmooth[Zsmooth <= 0] = np.nan

    return Zsmooth


def log_smooth_plot(ax, energy, angle, Z, title, clim):
    Zplot = smooth_for_plot(Z)

    im = ax.imshow(
        Zplot,
        origin="lower",
        aspect="auto",
        extent=[energy.min(), energy.max(), angle.min(), angle.max()],
        cmap="turbo",
        norm=LogNorm(vmin=clim[0], vmax=clim[1]),
    )

    ax.set_xlim(X_LIMIT)
    ax.set_ylim(Y_LIMIT)

    ax.set_xlabel("Incident Energy [eV]", fontsize=13, fontweight="bold")
    ax.set_ylabel("Incidence Angle [deg from normal]", fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=12)

    plt.colorbar(im, ax=ax, label="IEAD [log scale]")

    return im


# ============================
# READ DATA
# ============================
with Dataset(file, "r") as nc:
    surfEDist_raw = nc.variables["surfEDist"][:]
    surfRefl_raw = nc.variables["surfReflDist"][:]

print("Raw surfEDist shape:", surfEDist_raw.shape)
print("Raw surfReflDist shape:", surfRefl_raw.shape)

# Python NetCDF order:
# surfEDist_raw = [bins, species, surface, energy, angle]
# Convert to:
# surfEDist = [angle, energy, surface, species, bins]
surfEDist = np.transpose(surfEDist_raw, (4, 3, 2, 1, 0))

# surfRefl_raw = [species, surface, energy, angle]
# Convert to:
# surfReflDist = [angle, energy, surface, species]
surfReflDist = np.transpose(surfRefl_raw, (3, 2, 1, 0))

print("Converted surfEDist shape [A,E,S,species,bins]:", surfEDist.shape)
print("Converted surfReflDist shape [A,E,S,species]:", surfReflDist.shape)

nA, nE, nSurf, nSpecies, nBins = surfEDist.shape

# ============================
# INCIDENT IEADs
# ============================
bins_to_plot = [BIN_COMBINED] + BIN_TA
labels = ["Combined", "Ta1+", "Ta2+", "Ta3+", "Ta4+"]

inc_AE_list = []

for b in bins_to_plot:
    if b > nBins:
        raise ValueError(f"Requested bin {b}, but only {nBins} bins exist.")

    bidx = b - 1

    # sum over species -> [angle, energy, surface]
    AES = np.nansum(surfEDist[:, :, :, :, bidx], axis=3)

    AE, angle_deg, energy_eV = prep_AES_to_AE(AES)

    inc_AE_list.append(AE)

Eidx = np.where((energy_eV >= E_MIN) & (energy_eV <= E_MAX))[0]

AE_comb = inc_AE_list[0]

# ============================
# INCIDENT 2x3 LOG-SCALE PLOT
# ============================
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(
    f"Incident IEADs at MPEX target, tilt = {tilt}°",
    fontsize=16,
    fontweight="bold",
)

positions = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)]

for k, pos in enumerate(positions):
    ax = axes[pos]
    log_smooth_plot(
        ax,
        energy_eV[Eidx],
        angle_deg,
        inc_AE_list[k][:, Eidx],
        labels[k],
        CLIM_INC_LOG,
    )

# 1D spectra
ax = axes[1, 0]

for k, label in enumerate(labels):
    spec = np.nanmean(inc_AE_list[k][:, Eidx], axis=0)
    ax.semilogy(energy_eV[Eidx], spec, linewidth=2, label=label)

ax.set_xlim(X_LIMIT)
ax.set_xlabel("Incident Energy [eV]", fontsize=13)
ax.set_ylabel(r"$\langle IEAD \rangle_\theta$", fontsize=13)
ax.set_title("Angle-averaged IEADs", fontsize=14, fontweight="bold")
ax.grid(True, which="both")
ax.legend()

fig.tight_layout()
fig.savefig("Incident_IEADs_2x3_log_smooth.png", dpi=600)

# ============================
# SEPARATE COMBINED INCIDENT LOG-SCALE PLOT
# ============================
fig, ax = plt.subplots(figsize=(9, 7))

log_smooth_plot(
    ax,
    energy_eV[Eidx],
    angle_deg,
    AE_comb[:, Eidx],
    f"Combined Incident IEAD, tilt = {tilt}°",
    CLIM_INC_LOG,
)

fig.tight_layout()
fig.savefig("Combined_Incident_IEAD_log_smooth.png", dpi=600)

# ============================
# REFLECTED EAD LOG-SCALE PLOT
# ============================
AES_refl = np.nansum(surfReflDist, axis=3)

refl_AE, refl_angle_deg, refl_energy_eV = prep_AES_to_AE(AES_refl)

Eidx_refl = np.where((refl_energy_eV >= E_MIN) & (refl_energy_eV <= E_MAX))[0]

refl_plot = refl_AE[:, Eidx_refl]

total_refl = np.nansum(refl_plot)
if total_refl > 0:
    refl_plot = refl_plot / total_refl

refl_plot = smooth_for_plot(refl_plot)

fig, ax = plt.subplots(figsize=(9, 7))

im = ax.imshow(
    refl_plot,
    origin="lower",
    aspect="auto",
    extent=[
        refl_energy_eV[Eidx_refl].min(),
        refl_energy_eV[Eidx_refl].max(),
        refl_angle_deg.min(),
        refl_angle_deg.max(),
    ],
    cmap="turbo",
    norm=LogNorm(vmin=CLIM_REFL_LOG[0], vmax=CLIM_REFL_LOG[1]),
)

ax.set_xlim(X_LIMIT)
ax.set_ylim(Y_LIMIT)

ax.set_xlabel("Reflected Energy [eV]", fontsize=13, fontweight="bold")
ax.set_ylabel("Reflection Angle [deg]", fontsize=13, fontweight="bold")
ax.set_title(f"Reflected EAD, tilt = {tilt}°", fontsize=14, fontweight="bold")

plt.colorbar(im, ax=ax, label="Reflected EAD [log scale]")

fig.tight_layout()
fig.savefig("Reflected_EAD_log_smooth.png", dpi=600)

plt.show()

print("Log-scale plotting complete.")