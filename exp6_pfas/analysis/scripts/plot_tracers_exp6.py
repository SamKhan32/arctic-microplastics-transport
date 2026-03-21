"""
plot_tracers_exp6.py
Plot tr1, tr2, tr3 transport side-by-side for EXP6_PFAS.
Produces per-timestep PNGs and optionally a GIF.

Usage (from EXP6_PFAS root):
    python plot_tracers_exp6.py [--gif]

Outputs: figs/tracerN_NNNNN.png, figs/tr_all_NNNNN.png
         figs/tr3_transport.gif  (if --gif)
"""

import sys
import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- config ---
RUN_DIR    = "archives/03_21_2026_arctic_output_exp6"
FIG_DIR    = "analysis/figures/postrun"
GRID_FILE  = "../shared/original_data/ECCO-GRID_06.nc"
BATHY_FILE = "../shared/inputs/bathy_arctic_open.bin"
NX, NY     = 270, 270
NR         = 1
DT_SNAP    = 1       # output frequency in days (adjust to match diags)

TRACERS = [
    ("PTRACER01", "Microplastics (tr1)", "Blues"),
    ("PTRACER02", "Gyre Tracer (tr2)",   "Greens"),
    ("PTRACER03", "PFAS (tr3)",          "Reds"),
]

# PFAS source markers for tr3 panel
PFAS_SOURCES = [
    (69.0, -135.0, "Mackenzie"),
    (66.0, -169.0, "Bering"),
    (79.0,    0.0, "Fram"),
]

# --- helpers ---
def load_grid():
    import netCDF4 as nc
    ds = nc.Dataset(GRID_FILE)
    try:
        lon2d = ds.variables["XC"][:]
        lat2d = ds.variables["YC"][:]
    except KeyError:
        lon2d = ds.variables["lonc"][:]
        lat2d = ds.variables["latc"][:]
    ds.close()
    return lon2d, lat2d

def load_bathy():
    b = np.fromfile(BATHY_FILE, dtype=">f8").reshape(NY, NX)
    return b

def find_tracer_files(varname):
    """Return sorted list of MITgcm binary diag files for varname."""
    pattern = os.path.join(RUN_DIR, f"{varname}.*.data")
    files = sorted(glob.glob(pattern))
    return files

def read_diag(path):
    """Read a single MITgcm float32 big-endian diagnostic file."""
    data = np.fromfile(path, dtype=">f4").reshape(NR, NY, NX)
    return data[0]   # surface layer

def compute_vmax(files, pct=95, max_frames=20):
    """Estimate robust vmax from a sample of frames."""
    sample = files[::max(1, len(files)//max_frames)][:max_frames]
    vals = []
    for f in sample:
        d = read_diag(f)
        if d.max() > 0:
            vals.append(np.percentile(d[d > 0], pct))
    return float(np.median(vals)) if vals else 1.0

def make_norm(vmax):
    return mcolors.LogNorm(vmin=1e-6 * vmax, vmax=vmax)

# --- main ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gif", action="store_true")
    parser.add_argument("--tr3only", action="store_true",
                        help="Only produce combined panel PNGs (faster)")
    args = parser.parse_args()

    os.makedirs(FIG_DIR, exist_ok=True)

    print("Loading grid and bathy...")
    lon2d, lat2d = load_grid()
    bathy = load_bathy()
    land_mask = bathy >= 0.0

    proj = ccrs.NorthPolarStereo(central_longitude=0)
    extent = [-180, 180, 55, 90]

    # --- collect files per tracer ---
    all_files = []
    vmaxes    = []
    for (varname, label, cmap) in TRACERS:
        files = find_tracer_files(varname)
        all_files.append(files)
        if files:
            vm = compute_vmax(files)
            print(f"  {varname}: {len(files)} frames, vmax={vm:.3e}")
        else:
            vm = 1.0
            print(f"  {varname}: NO FILES FOUND (check RUN_DIR and diags config)")
        vmaxes.append(vm)

    n_frames = max(len(f) for f in all_files) if any(all_files) else 0
    if n_frames == 0:
        print("ERROR: no diagnostic files found in", RUN_DIR)
        sys.exit(1)

    gif_frames = []

    for frame_idx in range(n_frames):
        fig, axes = plt.subplots(
            1, 3, figsize=(18, 7),
            subplot_kw={"projection": proj}
        )
        fig.patch.set_facecolor("0.12")

        for ax_idx, (ax, (varname, label, cmap)) in enumerate(zip(axes, TRACERS)):
            files = all_files[ax_idx]
            norm  = make_norm(vmaxes[ax_idx])

            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND, facecolor="0.3", zorder=2)
            ax.coastlines(resolution="50m", color="0.7", linewidth=0.5, zorder=3)
            ax.gridlines(color="0.4", linewidth=0.3, linestyle="--")

            if frame_idx < len(files):
                data = read_diag(files[frame_idx])
                data = np.where(land_mask, np.nan, data)
                data = np.where(data <= 0, np.nan, data)

                pcm = ax.pcolormesh(
                    lon2d, lat2d, data,
                    norm=norm, cmap=cmap,
                    transform=ccrs.PlateCarree(),
                    shading="auto", zorder=1
                )
                cb = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
                cb.ax.tick_params(labelsize=7, colors="0.8")
                cb.set_label("Tracer conc.", fontsize=8, color="0.8")
            else:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", color="0.6")

            # mark PFAS source points on tr3 panel
            if ax_idx == 2:
                for (slat, slon, slabel) in PFAS_SOURCES:
                    ax.plot(slon, slat, "y*", ms=8,
                            transform=ccrs.PlateCarree(), zorder=5)
                    ax.text(slon + 3, slat, slabel, fontsize=6, color="yellow",
                            transform=ccrs.PlateCarree(), zorder=5)

            day = frame_idx * DT_SNAP
            ax.set_title(f"{label}\nDay {day}", fontsize=9, color="white", pad=4)

        plt.tight_layout(pad=1.2)
        out_path = os.path.join(FIG_DIR, f"tr_all_{frame_idx:05d}.png")
        fig.savefig(out_path, dpi=120, facecolor=fig.get_facecolor())
        plt.close(fig)

        if args.gif:
            gif_frames.append(out_path)

        if frame_idx % 10 == 0:
            print(f"  frame {frame_idx}/{n_frames}")

    print(f"Done. PNGs in {FIG_DIR}/")

    if args.gif and gif_frames:
        import imageio
        gif_path = os.path.join(FIG_DIR, "tr3_transport.gif")
        imgs = [imageio.imread(f) for f in gif_frames]
        imageio.mimsave(gif_path, imgs, fps=6, loop=0)
        print(f"GIF written: {gif_path}")

if __name__ == "__main__":
    main()