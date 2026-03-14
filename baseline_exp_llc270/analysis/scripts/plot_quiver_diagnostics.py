"""
plot_quiver_gif.py
------------------
Animated quiver GIF of Arctic circulation (speed background + flow arrows).
One frame per output timestep.

Run from experiment root:
    python analysis/scripts/plot_quiver_gif.py

Output: analysis/figures/arctic_quiver.gif
"""

import os
import re
import sys
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
from netCDF4 import Dataset

# -----------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------
DELTA_T        = 1200   # seconds per timestep
STRIDE         = 3      # plot every Nth vector
SCALE          = 1.0    # quiver scale; lower = longer arrows
FRAME_DURATION = 0.5    # seconds per GIF frame

# -----------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
OUTPUT_DIR  = os.path.join(PROJECT_ROOT, "output")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "analysis", "figures")
ECCO_NC     = os.path.join(PROJECT_ROOT, "original_data", "ECCO-GRID_06.nc")

NX, NY = 270, 270
DTYPE  = ">f4"

# -----------------------------------------------------------------------

def read_field(field, timestep):
    fpath = os.path.join(OUTPUT_DIR, f"{field}.{timestep:010d}.data")
    arr = np.fromfile(fpath, dtype=DTYPE)
    return arr.reshape(NY, NX)

def list_timesteps(field="U"):
    pattern = re.compile(rf"^{field}\.(\d{{10}})\.data$")
    ts = [int(m.group(1)) for f in os.listdir(OUTPUT_DIR)
          if (m := pattern.match(f))]
    return sorted(t for t in ts if t > 0)

def mask_dry(arr):
    out = arr.copy().astype(float)
    out[arr == 0.0] = np.nan
    return out

# -----------------------------------------------------------------------

def render_frame(xc, yc, u, v, spd, ts, vmax):
    day = ts * DELTA_T / 86400.0

    sl = (slice(None, None, STRIDE), slice(None, None, STRIDE))
    xq, yq = xc[sl], yc[sl]
    uq, vq = u[sl], v[sl]

    proj     = ccrs.NorthPolarStereo(central_longitude=0)
    data_crs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(7, 7), dpi=120)
    ax  = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent([-180, 180, 60, 90], crs=data_crs)
    fig.patch.set_facecolor("#0d0d0d")
    ax.set_facecolor("#0d0d0d")

    pc = ax.pcolormesh(
        xc, yc, spd,
        cmap="plasma", vmin=0, vmax=vmax,
        shading="auto", transform=data_crs, rasterized=True
    )

    ax.quiver(
        xq, yq, uq, vq,
        scale=SCALE,
        scale_units="width",
        width=0.002,
        color="white",
        alpha=0.7,
        transform=data_crs,
        regrid_shape=20
    )
    land = cfeature.NaturalEarthFeature(
        'physical', 'land', '110m',
        facecolor='#2a2a2a', edgecolor='none')
    ax.add_feature(land, zorder=2)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, color="white", zorder=3)
    ax.gridlines(linewidth=0.4, color="white", alpha=0.3, linestyle="--")

    cbar = plt.colorbar(pc, ax=ax, fraction=0.035, pad=0.04, shrink=0.7)
    cbar.set_label("Speed (m/s)", color="white", fontsize=9)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    ax.set_title(
        f"Arctic Circulation  |  Day {day:.0f}",
        color="white", fontsize=11, pad=10
    )

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    frame = imageio.imread(buf)
    buf.close()
    return frame[:, :, :3]

# -----------------------------------------------------------------------

def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    timesteps = list_timesteps("U")
    if not timesteps:
        print(f"ERROR: no U files found in {OUTPUT_DIR}")
        sys.exit(1)
    print(f"Found {len(timesteps)} timesteps")

    print("Loading grid...")
    with Dataset(ECCO_NC) as ds:
        xc = np.array(ds.variables["XC"][:]).squeeze().reshape(NY, NX)
        yc = np.array(ds.variables["YC"][:]).squeeze().reshape(NY, NX)

    print("Computing global color scale from all timesteps...")
    all_spd = []
    for ts in timesteps:
        u = mask_dry(read_field("U", ts))
        v = mask_dry(read_field("V", ts))
        all_spd.append(np.nanpercentile(np.sqrt(u**2 + v**2), 95))
    vmax = float(np.median(all_spd))
    print(f"  vmax = {vmax:.4e} m/s")

    frames = []
    for i, ts in enumerate(timesteps):
        day = ts * DELTA_T / 86400.0
        u   = mask_dry(read_field("U", ts))
        v   = mask_dry(read_field("V", ts))
        spd = np.sqrt(u**2 + v**2)
        spd[spd == 0.0] = np.nan
        print(f"  Frame {i+1}/{len(timesteps)}  day {day:.0f}  "
              f"speed max={np.nanmax(spd):.3e} m/s")
        frames.append(render_frame(xc, yc, u, v, spd, ts, vmax))

    gif_path = os.path.join(FIGURES_DIR, "arctic_quiver.gif")
    print(f"Writing GIF -> {gif_path}")
    imageio.mimwrite(gif_path, frames, duration=FRAME_DURATION, loop=0)
    print(f"Done. {len(frames)} frames.")

if __name__ == "__main__":
    main()