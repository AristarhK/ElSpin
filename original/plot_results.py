"""Plot exported CSV data. Dependencies: numpy, matplotlib. No simulation is altered."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", default=".")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--out", default="results.png")
    args = parser.parse_args()
    directory = Path(args.directory)
    status = (directory / (args.prefix + "run_status.txt")).read_text()
    if not status.startswith("NUMERICALLY_CONVERGED"):
        raise RuntimeError("The last run did not finish; inspect run_status.txt.")
    jet = np.genfromtxt(directory / (args.prefix + "jet_profile.csv"), delimiter=",", names=True)
    field = np.genfromtxt(directory / (args.prefix + "potential_rz.csv"), delimiter=",", names=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.18})
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    z = jet["z_m"] * 1000
    datasets = [("diameter_m", 1e6, "Jet diameter, um", "#176f9e"),
                ("V_m_per_s", 1, "Axial speed, m/s", "#146c59"),
                ("T_C", 1, "Temperature, deg C", "#bd5437"),
                ("polymer_mass_fraction", 100, "Polymer mass fraction, %", "#7054a5")]
    for ax, (key, scale, label, color) in zip(axes.ravel(), datasets):
        ax.plot(z, jet[key]*scale, color=color, lw=2)
        ax.set(xlabel="Distance from emitter, mm", ylabel=label, xscale="symlog")
        ax.set_xscale("symlog", linthresh=0.02)
    ax = axes[1, 1]
    ax.plot(z, jet["I_conv_A"]*1e9, label="Convective", lw=2, color="#176f9e")
    ax.plot(z, jet["I_cond_A"]*1e9, label="Ohmic", lw=2, color="#b27d18")
    ax.plot(z, jet["I_total_A"]*1e9, "--", label="Total", color="#333333")
    ax.set(xlabel="Distance from emitter, mm", ylabel="Current, nA")
    ax.set_xscale("symlog", linthresh=0.02)
    ax.legend(frameon=False)
    r, zz = np.unique(field["r_m"]), np.unique(field["z_m"])
    potential = field["Phi_V"].reshape(len(zz), len(r))
    ax = axes[1, 2]
    near_r, near_z = r <= 0.003, zz <= 0.003
    ax.grid(False)
    im = ax.pcolormesh(zz[near_z]*1e3, r[near_r]*1e3,
                       potential[np.ix_(near_z, near_r)].T / 1000, shading="auto", cmap="viridis")
    ax.set(xlabel="z, mm", ylabel="r, mm", title="Near-emitter potential")
    fig.colorbar(im, ax=ax, label="Potential, kV")
    fig.suptitle("Electrospinning | illustrative reduced-model calculation", fontsize=17, fontweight="bold")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, facecolor="white")
    print(out.resolve())


if __name__ == "__main__":
    main()
