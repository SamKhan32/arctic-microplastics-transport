"""
plot_ptracer_animation.py
Animate PTRACER01 and PTRACER02 overlaid on the LLC270 Arctic cap grid.
Run from the exp4_wind_fix experiment directory.

Output: analysis/figures/postrun/ptracer_animation.gif
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize, LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import os

# --- paths ---
MITGRID  = 'run/input/arctic_cap.mitgrid'
DATA_DIR = 'archive/arctic_output_exp4'
OUT_DIR  = 'analysis/figures/postrun'
OUT_FILE = os.path.join(OUT_DIR, 'ptracer_animation.gif')

NX, NY   = 270, 270
N_FIELDS = 16
DT       = 1200

# --- load grid ---
print('Loading grid...')
d      = np.fromfile(MITGRID, dtype='>f8')
fields = d.reshape(N_FIELDS, NY, NX)
XC     = fields[0].ravel()
YC     = fields[1].ravel()

# --- load snapshots ---
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

def timestep_from_file(f):
    base = os.path.basename(f)
    for p in base.split('.'):
        if p.isdigit() and len(p) == 10:
            return int(p)
    return 0

timesteps = [timestep_from_file(f) for f in tr1_files[:n_frames]]

# --- fixed color normalization across all frames ---
all_tr1 = np.concatenate([s for s in tr1_snaps])
all_tr2 = np.concatenate([s for s in tr2_snaps])

tr1_vmax = float(np.percentile(all_tr1[all_tr1 > 1e-8], 99)) if np.any(all_tr1 > 1e-8) else 1.0
tr2_vmax = float(np.percentile(all_tr2[all_tr2 > 1e-8], 99)) if np.any(all_tr2 > 1e-8) else 1.0
tr1_vmin = 1e-6
tr2_vmin = 1e-6

print(f'TR1 color range: {tr1_vmin:.2e} to {tr1_vmax:.2e}')
print(f'TR2 color range: {tr2_vmin:.2e} to {tr2_vmax:.2e}')

norm1 = LogNorm(vmin=tr1_vmin, vmax=tr1_vmax)
norm2 = LogNorm(vmin=tr2_vmin, vmax=tr2_vmax)

cmap1 = plt.cm.Blues
cmap2 = plt.cm.Oranges

# --- setup figure ---
proj     = ccrs.NorthPolarStereo(central_longitude=0)
data_crs = ccrs.PlateCarree()

fig = plt.figure(figsize=(8, 8), facecolor='#0a0a1a')
ax  = fig.add_subplot(1, 1, 1, projection=proj, facecolor='#0a0a1a')
ax.set_extent([-180, 180, 60, 90], crs=data_crs)
ax.add_feature(cfeature.OCEAN,     facecolor='#0d1b2a', zorder=0)
ax.add_feature(cfeature.LAND,      facecolor='#1a1a2e', zorder=1)
ax.add_feature(cfeature.COASTLINE, edgecolor='#444466', linewidth=0.5, zorder=2)
ax.gridlines(color='#333355', linewidth=0.4, linestyle='--', zorder=2)

# seed first frame so scatter objects are properly initialized
def make_offsets(mask):
    if mask.any():
        return np.column_stack([XC[mask], YC[mask]])
    return np.zeros((1, 2))

def make_colors(vals, mask, norm, cmap):
    if mask.any():
        return cmap(norm(vals[mask]))
    return np.array([[0, 0, 0, 0]])

tr1_0   = tr1_snaps[0]
tr2_0   = tr2_snaps[0]
mask1_0 = tr1_0 > tr1_vmin
mask2_0 = tr2_0 > tr2_vmin

sc1 = ax.scatter(
    XC[mask1_0] if mask1_0.any() else [0],
    YC[mask1_0] if mask1_0.any() else [0],
    s=4,
    c=norm1(tr1_0[mask1_0]) if mask1_0.any() else [0],
    cmap=cmap1, norm=norm1,
    alpha=0.85, linewidths=0,
    transform=data_crs, zorder=3, label='TR1 microplastics')

sc2 = ax.scatter(
    XC[mask2_0] if mask2_0.any() else [0],
    YC[mask2_0] if mask2_0.any() else [0],
    s=12,
    c=norm2(tr2_0[mask2_0]) if mask2_0.any() else [0],
    cmap=cmap2, norm=norm2,
    alpha=0.95, linewidths=0,
    transform=data_crs, zorder=4, label='TR2 gyre marker')

title = ax.set_title('', color='white', fontsize=11, pad=10)
ax.legend(loc='lower left', fontsize=8, framealpha=0.3,
          labelcolor='white', markerscale=3, facecolor='#0a0a1a')

cb1 = fig.colorbar(plt.cm.ScalarMappable(norm=norm1, cmap=cmap1),
                   ax=ax, fraction=0.025, pad=0.04, shrink=0.5)
cb1.set_label('TR1 concentration', color='white', fontsize=8)
cb1.ax.yaxis.set_tick_params(color='white')
plt.setp(cb1.ax.yaxis.get_ticklabels(), color='white')
cb1.outline.set_edgecolor('white')

cb2 = fig.colorbar(plt.cm.ScalarMappable(norm=norm2, cmap=cmap2),
                   ax=ax, fraction=0.025, pad=0.08, shrink=0.5)
cb2.set_label('TR2 concentration', color='white', fontsize=8)
cb2.ax.yaxis.set_tick_params(color='white')
plt.setp(cb2.ax.yaxis.get_ticklabels(), color='white')
cb2.outline.set_edgecolor('white')

fig.patch.set_facecolor('#0a0a1a')

def update(frame):
    tr1   = tr1_snaps[frame]
    tr2   = tr2_snaps[frame]
    mask1 = tr1 > tr1_vmin
    mask2 = tr2 > tr2_vmin

    if mask1.any():
        sc1.set_offsets(np.column_stack([XC[mask1], YC[mask1]]))
        sc1.set_array(tr1[mask1])
    else:
        sc1.set_offsets(np.zeros((1, 2)))
        sc1.set_array(np.array([tr1_vmin]))

    if mask2.any():
        sc2.set_offsets(np.column_stack([XC[mask2], YC[mask2]]))
        sc2.set_array(tr2[mask2])
    else:
        sc2.set_offsets(np.zeros((1, 2)))
        sc2.set_array(np.array([tr2_vmin]))

    days = timesteps[frame] * DT / 86400
    title.set_text(f'Arctic tracer transport  |  day {days:.0f}  (step {timesteps[frame]})')
    return sc1, sc2, title

ani = animation.FuncAnimation(fig, update, frames=n_frames,
                               interval=100, blit=False)

os.makedirs(OUT_DIR, exist_ok=True)
print(f'Saving to {OUT_FILE} ...')
ani.save(OUT_FILE, writer='pillow', fps=10, dpi=120)
print('Done.')
plt.close(fig)