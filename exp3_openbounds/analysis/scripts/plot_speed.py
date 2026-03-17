"""
plot_uv_speed.py
Animate current speed (sqrt(U^2 + V^2)) on the LLC270 Arctic cap grid.
Check for gradient structure near open boundaries over full run.
Run from the exp2_ptracers experiment directory.

Output: analysis/figures/postrun/uv_speed.gif
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import os

# --- paths ---
MITGRID  = 'run/input/arctic_cap.mitgrid'
DATA_DIR = 'archive/arctic_output_exp3_coeff_fixed'
OUT_DIR  = 'analysis/figures/postrun'
OUT_FILE = os.path.join(OUT_DIR, 'uv_speed.gif')

NX, NY   = 270, 270
N_FIELDS = 16
DT       = 1200  # seconds per timestep

# --- load grid ---
print('Loading grid...')
d = np.fromfile(MITGRID, dtype='>f8')
fields = d.reshape(N_FIELDS, NY, NX)
XC = fields[0].ravel()
YC = fields[1].ravel()

# --- load all snapshots ---
def load_snapshots(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError('No files matched: ' + pattern)
    snaps = []
    for f in files:
        snaps.append(np.fromfile(f, dtype='>f4').reshape(NY, NX))
    print(f'  loaded {len(snaps)} snapshots from {pattern}')
    return snaps, files

def timestep_from_file(f):
    base = os.path.basename(f)
    for p in base.split('.'):
        if p.isdigit() and len(p) == 10:
            return int(p)
    return 0

print('Loading U...')
u_snaps, u_files = load_snapshots(os.path.join(DATA_DIR, 'U.*.data'))
print('Loading V...')
v_snaps, _ = load_snapshots(os.path.join(DATA_DIR, 'V.*.data'))

n_frames  = min(len(u_snaps), len(v_snaps))
timesteps = [timestep_from_file(f) for f in u_files[:n_frames]]
print(f'Animating {n_frames} frames...')

# compute speed for all frames, get global vmax from 99th percentile
all_speeds = [np.sqrt(u_snaps[i]**2 + v_snaps[i]**2).ravel() for i in range(n_frames)]
nonzero    = np.concatenate([s[s > 0] for s in all_speeds])
vmax       = float(np.percentile(nonzero, 99)) if len(nonzero) > 0 else 1e-4

# --- plot ---
proj     = ccrs.NorthPolarStereo(central_longitude=0)
data_crs = ccrs.PlateCarree()

os.makedirs(OUT_DIR, exist_ok=True)

fig = plt.figure(figsize=(8, 8), facecolor='#0a0a1a')
ax  = fig.add_subplot(1, 1, 1, projection=proj, facecolor='#0a0a1a')
ax.set_extent([-180, 180, 60, 90], crs=data_crs)

ax.add_feature(cfeature.OCEAN,     facecolor='#0d1b2a', zorder=0)
ax.add_feature(cfeature.LAND,      facecolor='#1a1a2e', zorder=1)
ax.add_feature(cfeature.COASTLINE, edgecolor='#444466', linewidth=0.5, zorder=2)
ax.gridlines(color='#333355', linewidth=0.4, linestyle='--', zorder=2)

sc = ax.scatter([], [], s=1, c=[], cmap='plasma',
                vmin=0, vmax=vmax, alpha=0.9, linewidths=0,
                transform=data_crs, zorder=3)

cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.04, shrink=0.6)
cb.set_label('speed (m/s)', color='white', fontsize=9)
cb.ax.yaxis.set_tick_params(color='white')
plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
cb.outline.set_edgecolor('white')

title = ax.set_title('', color='white', fontsize=11, pad=10)
fig.patch.set_facecolor('#0a0a1a')

def update(frame):
    speed = all_speeds[frame]
    sc.set_offsets(np.column_stack([XC, YC]))
    sc.set_array(speed)
    days = timesteps[frame] * DT / 86400
    title.set_text(f'Current speed  |  day {days:.0f}  (step {timesteps[frame]})')
    return sc, title

ani = animation.FuncAnimation(fig, update, frames=n_frames,
                               interval=100, blit=False)

print(f'Saving to {OUT_FILE} ...')
ani.save(OUT_FILE, writer='pillow', fps=10, dpi=120)
print('Done.')
plt.close(fig)