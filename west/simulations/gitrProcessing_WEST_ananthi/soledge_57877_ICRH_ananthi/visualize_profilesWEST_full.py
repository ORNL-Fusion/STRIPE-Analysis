#!/usr/bin/env python3
"""
Full visualization script for profilesWEST.nc.

Reads:
  - profilesWEST.nc

Optionally reads:
  - ICRH_57877_8s.mat for experimental OMP reflectometry comparison

Plots:
  - B fields: br, bt, bz
  - Electrons: ne, te
  - D+: ni, ti, vr, vt, vz, vp
  - Oxygen O1+--O8+: no1--no8, to1--to8
  - Optional oxygen velocities: vopar_o#, vro#, vto#, vzo#, vpo#
  - Geometry stored in NetCDF:
      wall_r, wall_z
      sep1_r, sep1_z
      sep2_r, sep2_z
  - Mid-Z and mid-R lineouts for ne and Te
  - Optional experimental vs simulation OMP ne comparison

This script does not modify profilesWEST.nc.
It does not save plots by default.
Antenna centroids are not overlaid.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Optional only needed for experimental comparison
try:
    from scipy.io import loadmat
    from scipy.io.matlab import mat_struct
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


# ============================================================
# User settings
# ============================================================

ncfile = "profilesWEST.nc"

outdir = Path("profilesWEST_visualization")
outdir.mkdir(exist_ok=True)

save_plots = False
show_plots = True
dpi = 300

# Plot switches
plot_b_fields = True
plot_electrons = True
plot_dplus = True
plot_oxygen_density = True
plot_oxygen_temperature = True
plot_oxygen_velocity = False  # creates many extra plots

plot_lineouts = True
compare_exp_omp_ne = True

# Oxygen charge states
nO = 8

# Experimental comparison settings
exp_ne_file = "ICRH_57877_8s.mat"
exp_t_min = 5.0
exp_t_max = 10.0
omp_z_target = 0.0

R_zoom = (2.9, 3.1)
ne_ylim = (1e16, 1e20)

# Density plot settings
use_log_density = True

# Matplotlib style
plt.rcParams.update({
    "font.size": 15,
    "axes.grid": True,
    "figure.dpi": 110,
})


# ============================================================
# General helpers
# ============================================================

def nc_var_exists(ds, name):
    return name in ds.variables


def read_nc_var(ds, name, required=False):
    if name not in ds.variables:
        if required:
            raise KeyError(f"Required NetCDF variable missing: {name}")
        print(f"Skipping missing variable: {name}")
        return None
    arr = np.asarray(ds.variables[name][:], dtype=float)
    return arr


def read_nc_vector(ds, name):
    arr = read_nc_var(ds, name, required=False)
    if arr is None:
        return np.array([])
    arr = np.ravel(arr).astype(float)
    arr = arr[np.isfinite(arr)]
    return arr


def orient_rz(A, R, Z, name="array"):
    """
    Ensure A has shape [nR, nZ].
    NetCDF usually stores [nR, nZ].
    """
    if A is None:
        return None

    A = np.asarray(A, dtype=float)
    A[~np.isfinite(A)] = 0.0

    nR = len(R)
    nZ = len(Z)

    if A.shape == (nR, nZ):
        return A
    if A.shape == (nZ, nR):
        return A.T

    raise ValueError(
        f"{name} has shape {A.shape}; expected {(nR, nZ)} or {(nZ, nR)}"
    )


def clean_positive(A):
    A = np.asarray(A, dtype=float).copy()
    A[~np.isfinite(A)] = np.nan
    A[A <= 0] = np.nan
    return A


def save_figure(fig, filename_base):
    if not save_plots:
        return
    png = outdir / f"{filename_base}.png"
    pdf = outdir / f"{filename_base}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"Saved: {png}")


# ============================================================
# MATLAB struct helpers for experimental comparison
# ============================================================

def _todict(obj):
    if not SCIPY_AVAILABLE:
        return obj
    if isinstance(obj, mat_struct):
        out = {}
        for name in obj._fieldnames:
            out[name] = _todict(getattr(obj, name))
        return out
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            return np.array([_todict(x) for x in obj], dtype=object)
        return obj
    return obj


def load_mat_struct(filename):
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for MATLAB .mat experimental comparison.")
    data = loadmat(filename, struct_as_record=False, squeeze_me=True)
    return {k: _todict(v) for k, v in data.items() if not k.startswith("__")}


def get_nested(dct, keys):
    val = dct
    for key in keys:
        val = val[key]
    return val


def ensure_2d(a):
    a = np.asarray(a, dtype=float)
    if a.ndim == 1:
        return a[:, None]
    return a


# ============================================================
# Read profilesWEST.nc
# ============================================================

if not Path(ncfile).is_file():
    raise FileNotFoundError(f"NetCDF file not found: {ncfile}")

with Dataset(ncfile, "r") as ds:
    R = read_nc_var(ds, "x", required=True).ravel()
    Z = read_nc_var(ds, "z", required=True).ravel()

    geom = {
        "wall_r": read_nc_vector(ds, "wall_r"),
        "wall_z": read_nc_vector(ds, "wall_z"),
        "sep1_r": read_nc_vector(ds, "sep1_r"),
        "sep1_z": read_nc_vector(ds, "sep1_z"),
        "sep2_r": read_nc_vector(ds, "sep2_r"),
        "sep2_z": read_nc_vector(ds, "sep2_z"),
    }

    variables = {}
    for name in [
        "br", "bt", "bz",
        "ne", "te",
        "ni", "ti",
        "vr", "vt", "vz", "vp",
    ]:
        arr = read_nc_var(ds, name, required=False)
        variables[name] = orient_rz(arr, R, Z, name) if arr is not None else None

    for q in range(1, nO + 1):
        for prefix in ["no", "to", "vopar_o", "vro", "vto", "vzo", "vpo"]:
            name = f"{prefix}{q}"
            arr = read_nc_var(ds, name, required=False)
            variables[name] = orient_rz(arr, R, Z, name) if arr is not None else None

print(f"Read {ncfile}: nR={len(R)}, nZ={len(Z)}")
print(f"R range: {np.nanmin(R):.6g} to {np.nanmax(R):.6g} m")
print(f"Z range: {np.nanmin(Z):.6g} to {np.nanmax(Z):.6g} m")


# ============================================================
# Plot helpers
# ============================================================

def overlay_geometry(ax):
    """Overlay wall and separatrix only. No antenna centroid overlay."""
    if geom["sep1_r"].size and geom["sep1_z"].size:
        ax.plot(geom["sep1_r"], geom["sep1_z"], "k--", linewidth=1.4)

    if geom["sep2_r"].size and geom["sep2_z"].size:
        ax.plot(geom["sep2_r"], geom["sep2_z"], "k--", linewidth=1.4)

    if geom["wall_r"].size and geom["wall_z"].size:
        ax.plot(geom["wall_r"], geom["wall_z"], "k-", linewidth=1.6)


def plot_2d(A, title, cbar_label, filename, log_scale=False):
    if A is None:
        return

    data = np.asarray(A, dtype=float)

    if log_scale:
        data_plot = clean_positive(data)
    else:
        data_plot = data.copy()
        data_plot[~np.isfinite(data_plot)] = np.nan

    fig, ax = plt.subplots(figsize=(8.8, 6.6))

    # pcolormesh wants [nZ, nR], so transpose [nR,nZ] -> [nZ,nR].
    mesh = ax.pcolormesh(R, Z, data_plot.T, shading="auto")

    if log_scale:
        vals = data_plot[np.isfinite(data_plot) & (data_plot > 0)]
        if vals.size:
            from matplotlib.colors import LogNorm
            mesh.set_norm(LogNorm(vmin=np.nanmin(vals), vmax=np.nanmax(vals)))

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)

    overlay_geometry(ax)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$R$ [m]")
    ax.set_ylabel(r"$Z$ [m]")
    ax.set_title(title)
    ax.set_xlim(np.nanmin(R), np.nanmax(R))
    ax.set_ylim(np.nanmin(Z), np.nanmax(Z))
    ax.grid(False)

    fig.tight_layout()
    save_figure(fig, filename)


# ============================================================
# 2D profile plots
# ============================================================

if plot_b_fields:
    plot_2d(variables["br"], r"$B_R$", r"$B_R$ [T]", "B_R", log_scale=False)
    plot_2d(variables["bt"], r"$B_t$", r"$B_t$ [T]", "B_t", log_scale=False)
    plot_2d(variables["bz"], r"$B_Z$", r"$B_Z$ [T]", "B_Z", log_scale=False)

if plot_electrons:
    plot_2d(variables["ne"], r"Electron density $n_e$", r"$n_e$ [m$^{-3}$]", "electron_density_ne", log_scale=use_log_density)
    plot_2d(variables["te"], r"Electron temperature $T_e$", r"$T_e$ [eV]", "electron_temperature_Te", log_scale=False)

if plot_dplus:
    plot_2d(variables["ni"], r"D$^+$ density $n_i$", r"$n_i$ [m$^{-3}$]", "Dplus_density_ni", log_scale=use_log_density)
    plot_2d(variables["ti"], r"D$^+$ temperature $T_i$", r"$T_i$ [eV]", "Dplus_temperature_Ti", log_scale=False)
    plot_2d(variables["vr"], r"D$^+$ velocity $V_R$", r"$V_R$ [m/s]", "Dplus_velocity_VR", log_scale=False)
    plot_2d(variables["vt"], r"D$^+$ velocity $V_t$", r"$V_t$ [m/s]", "Dplus_velocity_Vt", log_scale=False)
    plot_2d(variables["vz"], r"D$^+$ velocity $V_Z$", r"$V_Z$ [m/s]", "Dplus_velocity_VZ", log_scale=False)
    plot_2d(variables["vp"], r"D$^+$ speed $|V|$", r"$|V|$ [m/s]", "Dplus_speed_Vp", log_scale=False)

for q in range(1, nO + 1):
    if plot_oxygen_density:
        plot_2d(
            variables.get(f"no{q}"),
            rf"O$^{{{q}+}}$ density",
            rf"$n_{{O^{{{q}+}}}}$ [m$^{{-3}}$]",
            f"O{q}plus_density",
            log_scale=use_log_density,
        )

    if plot_oxygen_temperature:
        plot_2d(
            variables.get(f"to{q}"),
            rf"O$^{{{q}+}}$ temperature",
            rf"$T_{{O^{{{q}+}}}}$ [eV]",
            f"O{q}plus_temperature",
            log_scale=False,
        )

    if plot_oxygen_velocity:
        plot_2d(variables.get(f"vopar_o{q}"), rf"O$^{{{q}+}}$ $V_\parallel$", r"$V_\parallel$ [m/s]", f"O{q}plus_vparallel", log_scale=False)
        plot_2d(variables.get(f"vro{q}"), rf"O$^{{{q}+}}$ $V_R$", r"$V_R$ [m/s]", f"O{q}plus_VR", log_scale=False)
        plot_2d(variables.get(f"vto{q}"), rf"O$^{{{q}+}}$ $V_t$", r"$V_t$ [m/s]", f"O{q}plus_Vt", log_scale=False)
        plot_2d(variables.get(f"vzo{q}"), rf"O$^{{{q}+}}$ $V_Z$", r"$V_Z$ [m/s]", f"O{q}plus_VZ", log_scale=False)
        plot_2d(variables.get(f"vpo{q}"), rf"O$^{{{q}+}}$ $|V|$", r"$|V|$ [m/s]", f"O{q}plus_Vp", log_scale=False)


# ============================================================
# Lineouts for ne and Te
# ============================================================

def plot_ne_te_lineouts():
    ne = variables.get("ne")
    te = variables.get("te")
    if ne is None or te is None:
        return

    z_mid_value = 0.5 * (np.nanmin(Z) + np.nanmax(Z))
    r_mid_value = 0.5 * (np.nanmin(R) + np.nanmax(R))

    iz_mid = int(np.nanargmin(np.abs(Z - z_mid_value)))
    ir_mid = int(np.nanargmin(np.abs(R - r_mid_value)))

    print(f"Lineout ne/Te vs R at Z={Z[iz_mid]:.6g} m, index={iz_mid}")
    print(f"Lineout ne/Te vs Z at R={R[ir_mid]:.6g} m, index={ir_mid}")

    fig, ax1 = plt.subplots(figsize=(8.8, 5.8))
    ax2 = ax1.twinx()

    ax1.plot(R, ne[:, iz_mid], linewidth=2.2, label=r"$n_e$")
    ax2.plot(R, te[:, iz_mid], linewidth=2.2, linestyle="--", label=r"$T_e$")

    ax1.set_xlabel(r"$R$ [m]")
    ax1.set_ylabel(r"$n_e$ [m$^{-3}$]")
    ax2.set_ylabel(r"$T_e$ [eV]")
    ax1.set_title(rf"Mid-Z lineout at $Z={Z[iz_mid]:.4f}$ m")
    ax1.grid(True)

    fig.tight_layout()
    save_figure(fig, "lineout_midZ_ne_Te_vs_R")

    fig, ax1 = plt.subplots(figsize=(8.8, 5.8))
    ax2 = ax1.twinx()

    ax1.plot(Z, ne[ir_mid, :], linewidth=2.2, label=r"$n_e$")
    ax2.plot(Z, te[ir_mid, :], linewidth=2.2, linestyle="--", label=r"$T_e$")

    ax1.set_xlabel(r"$Z$ [m]")
    ax1.set_ylabel(r"$n_e$ [m$^{-3}$]")
    ax2.set_ylabel(r"$T_e$ [eV]")
    ax1.set_title(rf"Mid-R lineout at $R={R[ir_mid]:.4f}$ m")
    ax1.grid(True)

    fig.tight_layout()
    save_figure(fig, "lineout_midR_ne_Te_vs_Z")


if plot_lineouts:
    plot_ne_te_lineouts()


# ============================================================
# Experimental OMP comparison
# ============================================================

def plot_expt_vs_sim_omp_ne():
    if not SCIPY_AVAILABLE:
        print("Skipping experimental comparison: scipy is not available.")
        return

    if not Path(exp_ne_file).is_file():
        print(f"Skipping experimental comparison: file not found: {exp_ne_file}")
        return

    ne = variables.get("ne")
    if ne is None:
        print("Skipping experimental comparison: ne not found in NetCDF.")
        return

    mat = load_mat_struct(exp_ne_file)
    reflec = get_nested(mat, ["data", "WDP", "S57877", "reflec"])

    t_exp = np.asarray(reflec["t"], dtype=float).ravel()
    idx_t = np.where((t_exp >= exp_t_min) & (t_exp <= exp_t_max))[0]

    if idx_t.size == 0:
        print(f"No experimental time points found in {exp_t_min}--{exp_t_max} s.")
        return

    r_exp_all = np.asarray(reflec["position"]["r"], dtype=float)
    ne_exp_all = np.asarray(reflec["ne"], dtype=float)

    r_exp_all = ensure_2d(r_exp_all)
    ne_exp_all = ensure_2d(ne_exp_all)

    # align columns with time
    if r_exp_all.shape[1] != t_exp.size and r_exp_all.shape[0] == t_exp.size:
        r_exp_all = r_exp_all.T
    if ne_exp_all.shape[1] != t_exp.size and ne_exp_all.shape[0] == t_exp.size:
        ne_exp_all = ne_exp_all.T

    if r_exp_all.shape[1] != t_exp.size:
        raise ValueError(f"Could not align r_exp with time: shape={r_exp_all.shape}, t={t_exp.size}")
    if ne_exp_all.shape[1] != t_exp.size:
        raise ValueError(f"Could not align ne_exp with time: shape={ne_exp_all.shape}, t={t_exp.size}")

    r_exp = r_exp_all[:, idx_t]
    ne_exp = ne_exp_all[:, idx_t]

    r_exp = np.asarray(r_exp, dtype=float)
    ne_exp = clean_positive(ne_exp)
    r_exp[~np.isfinite(r_exp)] = np.nan

    print(
        f"Experimental profiles used: {idx_t.size} time slices, "
        f"t={np.nanmin(t_exp[idx_t]):.3f}--{np.nanmax(t_exp[idx_t]):.3f} s"
    )

    r_mean = np.nanmean(r_exp, axis=1)
    ne_mean = np.nanmean(ne_exp, axis=1)
    ne_std = np.nanstd(ne_exp, axis=1)
    good_exp = np.isfinite(r_mean) & np.isfinite(ne_mean) & (ne_mean > 0)

    # Smooth experimental mean by interpolating each time slice onto common R
    Rmin_exp = np.nanmin(r_exp)
    Rmax_exp = np.nanmax(r_exp)
    R_smooth = np.linspace(Rmin_exp, Rmax_exp, 300)

    ne_interp = np.full((R_smooth.size, idx_t.size), np.nan)

    for j in range(idx_t.size):
        rr = r_exp[:, j]
        nn = ne_exp[:, j]
        good = np.isfinite(rr) & np.isfinite(nn) & (nn > 0)

        if np.count_nonzero(good) >= 2:
            rr_good = rr[good]
            nn_good = nn[good]

            order = np.argsort(rr_good)
            rr_good = rr_good[order]
            nn_good = nn_good[order]

            rr_unique, unique_idx = np.unique(rr_good, return_index=True)
            nn_unique = nn_good[unique_idx]

            if rr_unique.size >= 2:
                ne_interp[:, j] = np.interp(
                    R_smooth,
                    rr_unique,
                    nn_unique,
                    left=np.nan,
                    right=np.nan,
                )

    ne_smooth_mean = np.nanmean(ne_interp, axis=1)
    ne_smooth_std = np.nanstd(ne_interp, axis=1)
    good_smooth = np.isfinite(R_smooth) & np.isfinite(ne_smooth_mean) & (ne_smooth_mean > 0)

    iz_omp = int(np.nanargmin(np.abs(Z - omp_z_target)))
    Z_omp = float(Z[iz_omp])
    ne_sim = ne[:, iz_omp]
    good_sim = np.isfinite(R) & np.isfinite(ne_sim) & (ne_sim > 0)

    print(f"Simulation OMP lineout: requested Z={omp_z_target:.4f} m, using Z={Z_omp:.6f} m, index={iz_omp}")

    fig, ax = plt.subplots(figsize=(9.5, 6.5))

    for j in range(ne_exp.shape[1]):
        ax.semilogy(r_exp[:, j], ne_exp[:, j], "-", color="0.80", linewidth=0.8)

    upper = ne_smooth_mean + ne_smooth_std
    lower = ne_smooth_mean - ne_smooth_std
    lower[lower <= 0] = np.nan

    good_band = (
        np.isfinite(R_smooth)
        & np.isfinite(upper)
        & np.isfinite(lower)
        & (upper > 0)
        & (lower > 0)
    )

    if np.any(good_band):
        ax.fill_between(
            R_smooth[good_band],
            lower[good_band],
            upper[good_band],
            alpha=0.35,
            label="Experiment ±1σ",
        )

    ax.errorbar(
        r_mean[good_exp],
        ne_mean[good_exp],
        yerr=ne_std[good_exp],
        fmt="ko",
        markersize=5,
        linewidth=1.4,
        markerfacecolor="k",
        label="Reflectometry average",
    )

    if np.any(good_smooth):
        ax.semilogy(
            R_smooth[good_smooth],
            ne_smooth_mean[good_smooth],
            linewidth=3.0,
            label="Reflectometry smooth average",
        )

    ax.semilogy(
        R[good_sim],
        ne_sim[good_sim],
        "r-",
        linewidth=3.0,
        label=f"SOLEDGE/GITR, Z={Z_omp:.4f} m",
    )

    if geom["wall_r"].size:
        ax.axvline(np.nanmin(geom["wall_r"]), color="k", linestyle=":", linewidth=1.2)
        ax.axvline(np.nanmax(geom["wall_r"]), color="k", linestyle=":", linewidth=1.2)

    ax.set_yscale("log")
    ax.set_xlim(*R_zoom)
    ax.set_ylim(*ne_ylim)
    ax.set_xlabel(r"$R$ [m]")
    ax.set_ylabel(r"$n_e$ [m$^{-3}$]")
    ax.set_title(f"OMP density comparison, {exp_t_min:.1f}--{exp_t_max:.1f} s")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend(loc="best")

    fig.tight_layout()
    save_figure(fig, "comparison_OMP_ne_expt_vs_sim_zoom_log")


if compare_exp_omp_ne:
    plot_expt_vs_sim_omp_ne()


if show_plots:
    plt.show()
else:
    plt.close("all")

print("Visualization complete.")
