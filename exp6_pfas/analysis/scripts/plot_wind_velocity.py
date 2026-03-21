# exp3_openbounds/analysis/scripts/plot_wind_velocity.py
# Two-panel diagnostic: rotated wind stress seasonal cycle + time-mean model velocity.
# Run from exp3_openbounds/

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob

# --- paths (relative to exp3_openbounds/) ---
GRID_NC    = '../shared/original_data/ECCO-GRID_06.nc'
TAUX_FILE  = '../shared/inputs/taux_monthly_rot.bin'
TAUY_FILE  = '../shared/inputs/tauy_monthly_rot.bin'
BATHY_FILE = '../shared/inputs/bathy_arctic_open.bin'
UVEL_DIR   = 'archive/arctic_output_exp4'
OUT_WIND   = 'analysis/figures/postrun/wind_stress_seasonal.png'
OUT_VEL    = 'analysis/figures/postrun/velocity_mean.png'

NX, NY = 270, 270
S = 10  # quiver skip
MONTHS = ['Jan','Feb','Mar','Apr','May','Jun',
          'Jul','Aug','Sep','Oct','Nov','Dec']

# --- load grid and mask ---
g       = nc.Dataset(GRID_NC)
XC      = np.array(g.variables['XC']).squeeze()
YC      = np.array(g.variables['YC']).squeeze()
bathy   = np.fromfile(BATHY_FILE, dtype='>f8').reshape(NY, NX)
ocean   = bathy < 0

# --- load wind stress ---
taux = np.fromfile(TAUX_FILE, dtype='>f8').reshape(12, NY, NX)
tauy = np.fromfile(TAUY_FILE, dtype='>f8').reshape(12, NY, NX)

# mask land
taux = np.where(ocean[np.newaxis], taux, np.nan)
tauy = np.where(ocean[np.newaxis], tauy, np.nan)

# --- figure 1: seasonal wind stress (3x4) ---
proj = ccrs.NorthPolarStereo()
fig, axes = plt.subplots(3, 4, subplot_kw={'projection': proj}, figsize=(18, 13))

# back-rotate grid-relative -> geographic for display
g = nc.Dataset(GRID_NC)
CS = np.array(g.variables['CS'][:]).squeeze()
SN = np.array(g.variables['SN'][:]).squeeze()
g.close()

for m, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180, 62, 90], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, color='lightgray', zorder=4)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, zorder=5)
    ax.gridlines(linewidth=0.2, color='gray', alpha=0.4)

    tx_rot = taux[m]
    ty_rot = tauy[m]
    tx_geo =  tx_rot * CS - ty_rot * SN
    ty_geo =  tx_rot * SN + ty_rot * CS

    speed = np.sqrt(tx_rot**2 + ty_rot**2)
    ax.pcolormesh(XC, YC, speed, transform=ccrs.PlateCarree(),
                  cmap='YlOrRd', shading='auto', zorder=2, vmin=0, vmax=0.15)
    ax.quiver(XC[::S, ::S], YC[::S, ::S],
              tx_geo[::S, ::S], ty_geo[::S, ::S],
              transform=ccrs.PlateCarree(),
              scale=0.8, width=0.004, color='black', zorder=3)
    ax.set_title(MONTHS[m], fontsize=9)

plt.suptitle('Rotated wind stress -- monthly climatology (ocean cells only)', fontsize=13)
plt.tight_layout()
plt.savefig(OUT_WIND, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_WIND}")

# --- load time-mean velocity ---
u_files = sorted(glob.glob(f'{UVEL_DIR}/U.*.data'))
v_files = sorted(glob.glob(f'{UVEL_DIR}/V.*.data'))
print(f"Loading {len(u_files)} velocity snapshots...")

u_sum = np.zeros((NY, NX))
v_sum = np.zeros((NY, NX))
for uf, vf in zip(u_files, v_files):
    u_sum += np.fromfile(uf, dtype='>f4').reshape(NY, NX)
    v_sum += np.fromfile(vf, dtype='>f4').reshape(NY, NX)
u_mean = np.where(ocean, u_sum / len(u_files), np.nan)
v_mean = np.where(ocean, v_sum / len(v_files), np.nan)
print(f"u_mean max={np.nanmax(np.abs(u_mean)):.4f}  v_mean max={np.nanmax(np.abs(v_mean)):.4f}")

# --- figure 2: time-mean velocity ---
fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(9, 9))
ax.set_extent([-180, 180, 62, 90], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, color='lightgray', zorder=4)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=5)
ax.gridlines(linewidth=0.3, color='gray', alpha=0.5)

speed = np.sqrt(u_mean**2 + v_mean**2)
ax.pcolormesh(XC, YC, speed, transform=ccrs.PlateCarree(),
              cmap='Blues', shading='auto', zorder=2)
ax.quiver(XC[::S, ::S], YC[::S, ::S],
          u_mean[::S, ::S], v_mean[::S, ::S],
          transform=ccrs.PlateCarree(),
          scale=0.3, width=0.004, color='black', zorder=3)
ax.set_title('Time-mean velocity -- 4-year mean (exp3)', fontsize=12)
plt.tight_layout()
plt.savefig(OUT_VEL, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_VEL}")