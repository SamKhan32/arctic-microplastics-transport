"""
plot_ptracer_animation.py
Animate PTRACER01 and PTRACER02 overlaid on the LLC270 Arctic cap grid.
Run from the exp2_ptracers experiment directory.

Output: analysis/figures/postrun/ptracer_animation.gif
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import os

# --- paths ---
MITGRID  = 'run/input/arctic_cap.mitgrid'
DATA_DIR = 'archive/2026_03_14_arctic_output_77'
OUT_DIR  = 'analysis/figures/postrun'
OUT_FILE = os.path.join(OUT_DIR, 'ptracer_animation.gif')

NX, NY   = 270, 270
N_FIELDS = 16
DT       = 1200  # seconds per timestep

# --- load grid coordinates ---
print('Loading grid...')
d = np.fromfile(MITGRID, dtype='>f8')
fields = d.reshape(N_FIELDS, NY, NX)
XC = fields[0].ravel()
YC = fields[1].ravel()

# --- load tracer snapshots ---
def load_snapshots(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError('No files matched: ' + pattern)
    snaps = []
    for f in files:
        arr = np.fromfile(f, dtype='>f4').reshape(NY, NX).ravel()
        snaps.append(arr)
    print(f'  loaded {len(snaps)} snapshots from {pattern}')
    return snaps, files

print('Loading PTRACER01...')
tr1_snaps, tr1_files = load_snapshots(os.path.join(DATA_DIR, 'PTRACER01.*.data'))
print('Loading PTRACER02...')
tr2_snaps, tr2_files = load_snapshots(os.path.join(DATA_DIR, 'PTRACER02.*.data'))

n_frames = min(len(tr1_snaps), len(tr2_snaps))
print(f'Animating {n_frames} frames...')

# --- timestep numbers ---
def timestep_from_file(f):
    base = os.path.basename(f)
    for p in base.split('.'):
        if p.isdigit() and len(p) == 10:
            return int(p)
    return 0

timesteps = [timestep_from_file(f) for f in tr1_files[:n_frames]]

# --- normalization ---
tr1_max = max(float(np.percentile(s[s > 0], 99)) for s in tr1_snaps if np.any(s > 0))
tr2_max = max(float(np.percentile(s[s > 0], 99)) for s in tr2_snaps if np.any(s > 0))

# --- projection ---
proj = ccrs.NorthPolarStereo(central_longitude=0)
data_crs = ccrs.PlateCarree()

os.makedirs(OUT_DIR, exist_ok=True)

fig = plt.figure(figsize=(8, 8), facecolor='#0a0a1a')
ax = fig.add_subplot(1, 1, 1, projection=proj, facecolor='#0a0a1a')

# set map extent to Arctic
ax.set_extent([-180, 180, 60, 90], crs=data_crs)

# map features
ax.add_feature(cfeature.OCEAN, facecolor='#0d1b2a', zorder=0)
ax.add_feature(cfeature.LAND, facecolor='#1a1a2e', zorder=1)
ax.add_feature(cfeature.COASTLINE, edgecolor='#444466', linewidth=0.5, zorder=2)
ax.gridlines(color='#333355', linewidth=0.4, linestyle='--', zorder=2)

# initial scatter plots
sc1 = ax.scatter([], [], s=2, c=[], cmap='Blues',
                 norm=Normalize(vmin=0, vmax=tr1_max),
                 alpha=0.85, linewidths=0,
                 transform=data_crs, zorder=3, label='TR1 microplastics')
sc2 = ax.scatter([], [], s=6, c=[], cmap='Oranges',
                 norm=Normalize(vmin=0, vmax=tr2_max),
                 alpha=0.95, linewidths=0,
                 transform=data_crs, zorder=4, label='TR2 gyre markers')

title = ax.set_title('', color='white', fontsize=11, pad=10)
ax.legend(loc='lower left', fontsize=8, framealpha=0.3,
          labelcolor='white', markerscale=5,
          facecolor='#0a0a1a')

cb1 = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, tr1_max), cmap='Blues'),
                   ax=ax, fraction=0.03, pad=0.04, shrink=0.6)
cb1.set_label('TR1 concentration', color='white', fontsize=8)
cb1.ax.yaxis.set_tick_params(color='white')
plt.setp(cb1.ax.yaxis.get_ticklabels(), color='white')
cb1.outline.set_edgecolor('white')

fig.patch.set_facecolor('#0a0a1a')

def update(frame):
    tr1 = tr1_snaps[frame]
    tr2 = tr2_snaps[frame]

    mask1 = tr1 > 1e-6
    mask2 = tr2 > 1e-6

    if mask1.any():
        sc1.set_offsets(np.column_stack([XC[mask1], YC[mask1]]))
        sc1.set_array(tr1[mask1])
    else:
        sc1.set_offsets(np.empty((0, 2)))

    if mask2.any():
        sc2.set_offsets(np.column_stack([XC[mask2], YC[mask2]]))
        sc2.set_array(tr2[mask2])
    else:
        sc2.set_offsets(np.empty((0, 2)))

    days = timesteps[frame] * DT / 86400
    title.set_text(f'Arctic tracer transport  |  day {days:.0f}  (step {timesteps[frame]})')
    return sc1, sc2, title

ani = animation.FuncAnimation(fig, update, frames=n_frames,
                               interval=100, blit=False)

print(f'Saving to {OUT_FILE} ...')
ani.save(OUT_FILE, writer='pillow', fps=10, dpi=120)
print('Done.')
plt.close(fig)