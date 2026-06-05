from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from utilities import compute_streamtube_stats


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
STREAMTUBE_FILE = "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462/streamtube_library.npz"
VELOCITY_FILE = "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/t092440_n014462/Run09_budget0_term01_t092440_n014462.s3D"

NX = 1389
NY = 631
NZ = 88

DX = 40/126
DY = 20/126
DZ = 10/126

X0 = 100.31746031563637
Y0 = 0.0
Z0 = 0.03968253968125

# Edit this directly, or replace it with np.loadtxt("x_stations.txt").
X_STATIONS = np.linspace(103, 300.0, 51)

CUTCELL_SAMPLES = (3, 3)
PERIODIC_X = False

OUTPUT_DATA = "S800_streamtube_stats.csv"
OUTPUT_FIGURE = "S800_streamtube_stats.png"


def write_stats(output, stats):
    output = Path(output)

    if output.suffix == ".npz":
        np.savez(output, **stats)
        return

    data = np.column_stack((stats["x"], stats["area"], stats["mass_flow_rate"]))
    np.savetxt(
        output,
        data,
        delimiter=",",
        header="x,area,mass_flow_rate",
        comments="",
    )


def plot_stats(stats, output):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True)

    axes[0].plot(stats["x"], stats["area"], color="tab:blue", linewidth=1.6)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("Area")
    axes[0].grid(True, alpha=0.25, linewidth=0.6)

    axes[1].plot(stats["x"], stats["mass_flow_rate"], color="tab:red", linewidth=1.6)
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("Mass flow rate")
    axes[1].grid(True, alpha=0.25, linewidth=0.6)

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    stats = compute_streamtube_stats(
        STREAMTUBE_FILE,
        VELOCITY_FILE,
        NX,
        NY,
        NZ,
        X_STATIONS,
        DX,
        DY,
        DZ,
        x0=X0,
        y0=Y0,
        z0=Z0,
        cutcell_samples=CUTCELL_SAMPLES,
        periodic_x=PERIODIC_X,
    )

    write_stats(OUTPUT_DATA, stats)
    plot_stats(stats, OUTPUT_FIGURE)

    print(OUTPUT_DATA)
    print(OUTPUT_FIGURE)


if __name__ == "__main__":
    main()
