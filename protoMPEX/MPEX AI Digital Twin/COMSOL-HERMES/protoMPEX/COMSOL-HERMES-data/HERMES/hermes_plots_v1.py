import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc


def make_dir(path):
    os.makedirs(path, exist_ok=True)


def save_png(fig, outdir, name):
    fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_case_name(filename):
    if "window-limited" in filename:
        return "window_limited"
    elif "limiter-limited" in filename:
        return "limiter_limited"
    else:
        return os.path.basename(os.path.dirname(filename)).replace("-", "_")


def plot_case(filename):
    case_name = get_case_name(filename)
    outdir = os.path.join("plots", case_name)
    make_dir(outdir)

    with nc.Dataset(filename, mode="r") as f:
        r = f.variables["r"][:]
        z = f.variables["z"][:]
        angle = f.variables["angle"][:]
        Ne = f.variables["electron_density"][:]
        Te = f.variables["electron_temperature"][:]

    X = r * np.cos(angle)
    Y = r * np.sin(angle)
    Z = z

    # -----------------------------
    # Save x, y, z, Ne, Te data
    # -----------------------------
    data = np.column_stack([
        X.ravel(),
        Y.ravel(),
        Z.ravel(),
        Ne.ravel(),
        Te.ravel()
    ])

    csv_file = os.path.join(outdir, f"{case_name}_xyz_ne_te.csv")
    np.savetxt(
        csv_file,
        data,
        delimiter=",",
        header="x[m],y[m],z[m],Ne[m^-3],Te[eV]",
        comments=""
    )

    print(f"Wrote data: {csv_file}")

    # -----------------------------
    # Z-slice X-Y plots
    # -----------------------------
    yind = min(20, Z.shape[1] - 1)

    fig, ax = plt.subplots()
    cm = ax.contourf(X[:, yind, :], Y[:, yind, :], Ne[:, yind, :], 50)
    plt.colorbar(cm, ax=ax, label="Ne [m$^{-3}$]")
    ax.set_aspect("equal")
    ax.set_title(f"{case_name}: density at Z = {Z[0, yind, 0]:.3f} m")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    save_png(fig, outdir, "density_xy_zslice")

    fig, ax = plt.subplots()
    cm = ax.contourf(X[:, yind, :], Y[:, yind, :], Te[:, yind, :], 50)
    plt.colorbar(cm, ax=ax, label="Te [eV]")
    ax.set_aspect("equal")
    ax.set_title(f"{case_name}: temperature at Z = {Z[0, yind, 0]:.3f} m")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    save_png(fig, outdir, "temperature_xy_zslice")

    # -----------------------------
    # theta = 0 r-z plots
    # -----------------------------
    theta_ind = np.unravel_index(np.argmin(np.abs(angle)), angle.shape)[2]

    R0 = r[:, :, theta_ind]
    Z0 = z[:, :, theta_ind]
    Ne0 = Ne[:, :, theta_ind]
    Te0 = Te[:, :, theta_ind]

    fig, ax = plt.subplots()
    cm = ax.contourf(R0, Z0, Ne0, 50)
    plt.colorbar(cm, ax=ax, label="Ne [m$^{-3}$]")
    ax.set_title(rf"{case_name}: density at $\theta = 0$")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    save_png(fig, outdir, "density_theta0_rz")

    fig, ax = plt.subplots()
    cm = ax.contourf(R0, Z0, Te0, 50)
    plt.colorbar(cm, ax=ax, label="Te [eV]")
    ax.set_title(rf"{case_name}: temperature at $\theta = 0$")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    save_png(fig, outdir, "temperature_theta0_rz")

    # -----------------------------
    # 1D r profiles at theta = 0
    # -----------------------------
    zind = Z0.shape[1] // 2

    fig, ax = plt.subplots()
    ax.plot(R0[:, zind], Ne0[:, zind], linewidth=2)
    ax.set_title(f"{case_name}: Ne radial profile at theta=0")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("Ne [m$^{-3}$]")
    ax.grid(True)
    save_png(fig, outdir, "ne_radial_profile_theta0")

    fig, ax = plt.subplots()
    ax.plot(R0[:, zind], Te0[:, zind], linewidth=2)
    ax.set_title(f"{case_name}: Te radial profile at theta=0")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("Te [eV]")
    ax.grid(True)
    save_png(fig, outdir, "te_radial_profile_theta0")

    # -----------------------------
    # 1D z profiles at theta = 0
    # -----------------------------
    rind = R0.shape[0] // 2

    fig, ax = plt.subplots()
    ax.plot(Z0[rind, :], Ne0[rind, :], linewidth=2)
    ax.set_title(f"{case_name}: Ne axial profile at theta=0")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("Ne [m$^{-3}$]")
    ax.grid(True)
    save_png(fig, outdir, "ne_axial_profile_theta0")

    fig, ax = plt.subplots()
    ax.plot(Z0[rind, :], Te0[rind, :], linewidth=2)
    ax.set_title(f"{case_name}: Te axial profile at theta=0")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("Te [eV]")
    ax.grid(True)
    save_png(fig, outdir, "te_axial_profile_theta0")

    # -----------------------------
    # 3D scatter plots
    # -----------------------------
    skip = max(1, X.size // 100000)

    Xs = X.ravel()[::skip]
    Ys = Y.ravel()[::skip]
    Zs = Z.ravel()[::skip]
    Nes = Ne.ravel()[::skip]
    Tes = Te.ravel()[::skip]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(Xs, Ys, Zs, c=Nes, s=2)
    fig.colorbar(p, ax=ax, label="Ne [m$^{-3}$]")
    ax.set_title(f"{case_name}: 3D density")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    save_png(fig, outdir, "density_3d")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(Xs, Ys, Zs, c=Tes, s=2)
    fig.colorbar(p, ax=ax, label="Te [eV]")
    ax.set_title(f"{case_name}: 3D temperature")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    save_png(fig, outdir, "temperature_3d")

    print(f"Finished case: {case_name}")
    print(f"Plots written to: {outdir}")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage:")
        print("  python hermes_plots_v1.py window-limited/260507/time_average.nc")
        print("  python hermes_plots_v1.py limiter-limited/260507/time_average.nc")
        sys.exit(1)

    filename = sys.argv[1]
    plot_case(filename)