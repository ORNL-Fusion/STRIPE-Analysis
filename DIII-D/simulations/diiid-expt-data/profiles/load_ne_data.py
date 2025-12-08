import pickle
import numpy as np
import matplotlib.pyplot as plt

def main():
    file_path = "profs_200882_1800_py3.pkl"

    # Load
    with open(file_path, "rb") as f:
        data = pickle.load(f)

    # Raw profile
    nedatpsi = data["nedatpsi"]
    psi_raw = np.array(nedatpsi["x"])
    ne_raw = np.array(nedatpsi["y"])
    ne_err = np.array(nedatpsi.get("err", np.ones_like(ne_raw) * 0.05))

    # Fitted profile
    netanhpsi = data["netanhpsi"]
    psi_fit = np.array(netanhpsi["x"])
    ne_fit = np.array(netanhpsi["y"])

    # Interpolate fitted profile onto raw psi grid
    ne_fit_on_raw = np.interp(psi_raw, psi_fit, ne_fit)

    # Print first few entries
    print("\nFirst 5 raw data points:")
    for i in range(5):
        print(f"psi={psi_raw[i]:.3f}, ne={ne_raw[i]:.3f}, err={ne_err[i]:.3f}, ne_fit_interp={ne_fit_on_raw[i]:.3f}")

    # Plot
    plt.figure(figsize=(7,5))

    # Raw data with error bars
    plt.errorbar(
        psi_raw,
        ne_raw,
        yerr=ne_err,
        fmt="o",
        capsize=3,
        label="Raw ne (nedatpsi)"
    )

    # Fitted profile evaluated at same psi
    plt.plot(
        psi_raw,
        ne_fit_on_raw,
        "r--",
        lw=2,
        label="Fitted ne (netanhpsi interpolated)"
    )

    plt.xlabel("Normalized poloidal flux $\\psi_N$")
    plt.ylabel("Electron density $n_e$ [$10^{20}$ m$^{-3}$]")
    plt.title("Electron Density Profile: Raw Data vs Fitted on Same Grid")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
