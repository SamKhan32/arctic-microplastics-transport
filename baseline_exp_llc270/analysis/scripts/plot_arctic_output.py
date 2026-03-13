"""
plot_arctic_output.py
---------------------
Produces a North Polar Stereographic animated GIF of current speed
from MITgcm single-CPU-IO binary output (LLC270 Arctic cap, barotropic).

Requires useSingleCpuIO=.TRUE. in PARM01 -- output files are:
    FIELD.0000000000.data  (single file per field per timestep)

Directory layout (run from experiment root):
    baseline_exp_llc270/
        output/              # .data/.meta files
        original_data/
            ECCO-GRID_06.nc
        analysis/
            scripts/         # this file
            figures/         # arctic_circulation.gif written here

Usage:
    python analysis/scripts/plot_arctic_output.py

Dependencies:
    numpy, matplotlib, cartopy, imageio, netCDF4
"""

import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
import io

# -----------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
OUTPUT_DIR  = os.path.join(PROJECT_ROOT, "output")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "analysis", "figures")
ECCO_NC     = os.path.join(PROJECT_ROOT, "original_data", "ECCO-GRID_06.nc")

# -----------------------------------------------------------------------
# GRID PARAMETERS
# -----------------------------------------------------------------------
NX    = 270
NY    = 270
DTYPE = ">f4"      # float32 big-endian (writeBinaryPrec default)

# GIF frame duration in seconds
FRAME_DURATION = 0.5

# -----------------------------------------------------------------------
# FILE I/O
# -----------------------------------------------------------------------

def read_field(output_dir, field, timestep):
    """
    Read a single-CPU-IO output file -> (NY, NX) float32 array.
    Filename format: FIELD.0000000000.data
    """
    fname = f"{field}.{timestep:010d}.data"
    fpath = os.path.join(output_dir, fname)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Missing file: {fname}\n  Looked in: {output_dir}")
    arr = np.fromfile(fpath, dtype=DTYPE)
    if arr.size != NX * NY:
        raise ValueError(
            f"{fname}: expected {NX*NY} values, got {arr.size}"
        )
    return arr.reshape(NY, NX)


def list_timesteps(output_dir, field="U"):
    """Return sorted list of integer timesteps available for a field."""
    pattern = re.compile(rf"^{field}\.(\d{{10}})\.data$")
    seen = set()
    for fname in os.listdir(output_dir):
        m = pattern.match(fname)
        if m:
            seen.add(int(m.group(1)))
    return sorted(seen)


def mask_dry(arr):
    """Replace exactly-zero cells with NaN (MITgcm dry cell fill value)."""
    out = arr.copy()
    out[arr == 0.0] = np.nan
    return out

# -----------------------------------------------------------------------
# GRID COORDINATES
# -----------------------------------------------------------------------

def load_grid_coords(ecco_nc):
    """
    Load XC, YC (tracer cell lon/lat) from ECCO-GRID_06.nc.
    Returns two (NY, NX) float arrays.
    """
    if not os.path.exists(ecco_nc):
        raise FileNotFoundError(f"ECCO grid file not found: {ecco_nc}")
    try:
        from netCDF4 import Dataset
        with Dataset(ecco_nc) as ds:
            xc = np.array(ds.variables["XC"][:]).squeeze()
            yc = np.array(ds.variables["YC"][:]).squeeze()
    except ImportError:
        import scipy.io.netcdf as snc
        with snc.netcdf_file(ecco_nc, "r") as ds:
            xc = ds.variables["XC"][:].copy().squeeze()
            yc = ds.variables["YC"][:].copy().squeeze()

    if xc.ndim == 2 and xc.shape != (NY, NX):
        tile_rows = xc.shape[0] // 13
        xc = xc[6 * tile_rows: 7 * tile_rows, :]
        yc = yc[6 * tile_rows: 7 * tile_rows, :]

    return xc.reshape(NY, NX), yc.reshape(NY, NX)

# -----------------------------------------------------------------------
# PLOTTING
# -----------------------------------------------------------------------

PROJECTION = ccrs.NorthPolarStereo(central_longitude=0)
DATA_CRS   = ccrs.PlateCarree()


def render_frame(xc, yc, speed, timestep, vmax, delta_t=300):
    """
    Render one GIF frame -> (H, W, 3) uint8 numpy array.
    """
    day = timestep * delta_t / 86400.0

    fig = plt.figure(figsize=(7, 7), dpi=120)
    ax = fig.add_subplot(1, 1, 1, projection=PROJECTION)
    ax.set_extent([-180, 180, 60, 90], crs=DATA_CRS)

    pc = ax.pcolormesh(
        xc, yc, speed,
        cmap="plasma",
        vmin=0, vmax=vmax,
        shading="auto",
        transform=DATA_CRS,
        rasterized=True
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, color="white", zorder=3)
    ax.add_feature(cfeature.LAND, facecolor="#2a2a2a", zorder=2)
    ax.gridlines(
        crs=DATA_CRS, draw_labels=False,
        linewidth=0.4, color="white", alpha=0.3, linestyle="--"
    )

    cbar = plt.colorbar(pc, ax=ax, fraction=0.035, pad=0.04, shrink=0.7)
    cbar.set_label("Current speed (m/s)", color="white", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    ax.set_title(
        f"Arctic Ocean Circulation\n"
        f"ERA5 January forcing  |  Day {day:.0f}",
        color="white", fontsize=11, pad=10
    )

    fig.patch.set_facecolor("#0d0d0d")
    ax.set_facecolor("#0d0d0d")

    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    frame = imageio.imread(buf)
    buf.close()
    return frame[:, :, :3]

# -----------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------

def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Output dir:   {OUTPUT_DIR}")

    print("Loading grid coordinates...")
    xc, yc = load_grid_coords(ECCO_NC)
    print(f"  Lon: {xc.min():.1f} to {xc.max():.1f}")
    print(f"  Lat: {yc.min():.1f} to {yc.max():.1f}")

    print("Scanning timesteps...")
    timesteps = list_timesteps(OUTPUT_DIR, field="U")
    if not timesteps:
        print(f"ERROR: no U files found in {OUTPUT_DIR}")
        sys.exit(1)
    timesteps = [ts for ts in timesteps if ts > 0]
    print(f"  Animating {len(timesteps)} frames: {timesteps}")

    print("Computing global color scale...")
    u0 = mask_dry(read_field(OUTPUT_DIR, "U", timesteps[0]))
    v0 = mask_dry(read_field(OUTPUT_DIR, "V", timesteps[0]))
    spd0 = np.sqrt(u0**2 + v0**2)
    vmax = float(np.nanpercentile(spd0, 99)) * 1.5
    print(f"  vmax = {vmax:.4e} m/s")

    frames = []
    for ts in timesteps:
        day = ts * 300 / 86400.0
        print(f"  Rendering Day {day:.0f} (timestep {ts})...")
        u   = mask_dry(read_field(OUTPUT_DIR, "U",   ts))
        v   = mask_dry(read_field(OUTPUT_DIR, "V",   ts))
        spd = np.sqrt(u**2 + v**2)
        spd[spd == 0.0] = np.nan
        print(
            f"    U max={np.nanmax(np.abs(u)):.3e}  "
            f"V max={np.nanmax(np.abs(v)):.3e}  "
            f"speed max={np.nanmax(spd):.3e} m/s"
        )
        frame = render_frame(xc, yc, spd, ts, vmax)
        frames.append(frame)

    gif_path = os.path.join(FIGURES_DIR, "arctic_circulation.gif")
    print(f"Writing GIF -> {gif_path}")
    imageio.mimwrite(gif_path, frames, duration=FRAME_DURATION, loop=0)
    print(f"Done. {len(frames)} frames at {FRAME_DURATION}s/frame.")


if __name__ == "__main__":
    main()