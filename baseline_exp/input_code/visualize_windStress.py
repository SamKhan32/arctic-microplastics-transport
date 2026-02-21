"""
visualize_wind_stress.py - Plot climatological monthly wind stress vectors
overlaid on bathymetry on an Arctic polar stereographic projection.

Produces one plot per month (or a selected month) showing:
  - Bathymetry as a filled contour background
  - Wind stress vectors subsampled onto a readable quiver grid
  - Stress magnitude as vector color

Usage:
    python visualize_wind_stress.py           # plots all 12 months
    python visualize_wind_stress.py --month 1 # January only
    python visualize_wind_stress.py --month 7 # July only
"""

import argparse
import os
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import make_grid, NX, NY

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT_DIR     = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MITGCM_DIR   = os.path.join(ROOT_DIR, "mitgcm_input")
ANALYSIS_DIR = os.path.join(ROOT_DIR, "analysis")
os.makedirs(ANALYSIS_DIR, exist_ok=True)

MONTH_NAMES = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
]

# Quiver subsampling: plot every Nth grid point in each direction.
# At NX=1530, NY=520, a stride of 30 gives ~50x17 arrows which is readable.
QUIVER_STRIDE = 10


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_bathy():
    path = os.path.join(MITGCM_DIR, "bathy.bin")
    if not os.path.exists(path):
        raise FileNotFoundError(f"bathy.bin not found at {path}")
    bathy = np.fromfile(path, dtype=">f4").reshape(NY, NX)
    return bathy


def load_stress(month):
    """Load tau_x, tau_y for a given month (1-indexed)."""
    fx = os.path.join(MITGCM_DIR, f"oceTauX_{month:02d}.bin")
    fy = os.path.join(MITGCM_DIR, f"oceTauY_{month:02d}.bin")
    for p in [fx, fy]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Not found: {p}. Run windStress.py first.")
    tau_x = np.fromfile(fx, dtype=">f4").reshape(NY, NX)
    tau_y = np.fromfile(fy, dtype=">f4").reshape(NY, NX)
    return tau_x, tau_y


# ---------------------------------------------------------------------------
# Plot one month
# ---------------------------------------------------------------------------
def plot_month(LON, LAT, bathy, month):
    tau_x, tau_y = load_stress(month)

    # Stress magnitude for coloring vectors
    magnitude = np.sqrt(tau_x * tau_x + tau_y * tau_y)

    # Mask land (bathy == 0) so it reads clearly
    bathy_masked = np.where(bathy == 0, np.nan, bathy)
    mag_masked   = np.where(bathy == 0, np.nan, magnitude)

    # Subsample for quiver
    sl = (slice(None, None, QUIVER_STRIDE), slice(None, None, QUIVER_STRIDE))
    qlon = LON[sl]
    qlat = LAT[sl]
    qtx  = tau_x[sl]
    qty  = tau_y[sl]
    qmag = mag_masked[sl]

    # Projection
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    data_crs = ccrs.PlateCarree()

    fig, ax = plt.subplots(
        figsize=(9, 9),
        subplot_kw={"projection": proj}
    )
    ax.set_extent([-180, 180, 55, 90], crs=data_crs)

    # Bathymetry background
    bathy_levels = np.linspace(-5500, 0, 56)
    cf = ax.contourf(
        LON, LAT, bathy_masked,
        levels=bathy_levels,
        cmap="Blues_r",
        transform=data_crs,
        extend="min"
    )

    # Land
    ax.add_feature(cfeature.LAND, facecolor="#c8b89a", zorder=3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="#5a4a3a", zorder=4)

    # Gridlines
    gl = ax.gridlines(
        draw_labels=False, linewidth=0.4,
        color="white", alpha=0.5, linestyle="--"
    )

    # Wind stress vectors, colored by magnitude
    vmax = 0.3  # N/m^2 — saturate colormap here; values above show as max color
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cmap_vec = plt.cm.plasma

    # Draw quiver — vectors are in grid-relative coords but we display them
    # in geographic space, so we need to transform back to east/north first.
    # We stored grid-relative (tau_x, tau_y); rotate back using the angle files.
    cos_a = np.fromfile(os.path.join(MITGCM_DIR, "angleCosC.bin"), dtype=">f4").reshape(NY, NX)
    sin_a = np.fromfile(os.path.join(MITGCM_DIR, "angleSinC.bin"), dtype=">f4").reshape(NY, NX)

    # Inverse rotation: tau_east = tau_x*cos - tau_y*sin
    #                   tau_north = tau_x*sin + tau_y*cos
    tau_e = tau_x * cos_a - tau_y * sin_a
    tau_n = tau_x * sin_a + tau_y * cos_a

    qte = tau_e[sl]
    qtn = tau_n[sl]

    # Ocean points only
    ocean = ~np.isnan(qmag)

    sc = ax.quiver(
        qlon[ocean], qlat[ocean],
        qte[ocean],  qtn[ocean],
        qmag[ocean],
        cmap=cmap_vec, norm=norm,
        transform=data_crs,
        scale=3.0,          # tune if arrows look too big or small
        scale_units="inches",
        width=0.003,
        headwidth=3,
        headlength=4,
        zorder=5
    )

    # Colorbars
    cbar_bathy = plt.colorbar(cf, ax=ax, fraction=0.025, pad=0.02, shrink=0.7)
    cbar_bathy.set_label("Depth (m)", fontsize=9)
    cbar_bathy.ax.tick_params(labelsize=8)

    cbar_vec = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap_vec),
        ax=ax, fraction=0.025, pad=0.06, shrink=0.7
    )
    cbar_vec.set_label("Stress magnitude (N/m^2)", fontsize=9)
    cbar_vec.ax.tick_params(labelsize=8)

    # Reference arrow
    ax.quiverkey(
        sc, X=0.12, Y=0.04, U=0.1,
        label="0.1 N/m^2", labelpos="E",
        fontproperties={"size": 8},
        coordinates="axes"
    )

    ax.set_title(
        f"ERA5 climatological wind stress - {MONTH_NAMES[month-1]}",
        fontsize=12, pad=10
    )

    out_path = os.path.join(ANALYSIS_DIR, f"wind_stress_{month:02d}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--month", type=int, default=None,
        help="Month to plot (1-12). Omit to plot all 12."
    )
    args = parser.parse_args()

    print("Loading grid ...")
    LON, LAT = make_grid()

    print("Loading bathymetry ...")
    bathy = load_bathy()

    months = [args.month] if args.month else list(range(1, 13))

    print(f"Plotting {len(months)} month(s) ...")
    for m in months:
        print(f"  {MONTH_NAMES[m-1]} ...", end=" ")
        plot_month(LON, LAT, bathy, m)

    print("Done. Plots written to analysis/")


if __name__ == "__main__":
    main()