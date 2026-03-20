# generate_wind.py
# Generates taux_monthly_rot.bin and tauy_monthly_rot.bin for MITgcm.
# Reads ERA5 monthly stress, computes 12-month climatology, interpolates
# onto LLC270 Arctic cap grid, rotates to grid-relative frame via CS/SN.
# Run from shared/

import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator
import os

ERA5_PATH = 'original_data/era5_stress_monthly_1991_2020.nc'
GRID_PATH = 'original_data/ECCO-GRID_06.nc'
OUT_DIR   = 'inputs'
os.makedirs(OUT_DIR, exist_ok=True)

# --- load ERA5 ---
ds  = nc.Dataset(ERA5_PATH)
lat = np.array(ds.variables['latitude'][:])   # descending: 90 to -90
lon = np.array(ds.variables['longitude'][:])  # 0 to 359.75
raw_u = np.array(ds.variables['avg_iews'][:]) # (360, 721, 1440) eastward stress
raw_v = np.array(ds.variables['avg_inss'][:]) # (360, 721, 1440) northward stress
ds.close()
print(f'ERA5 loaded: u shape={raw_u.shape}  lat[0]={lat[0]:.1f}')

# flip to ascending latitude for RegularGridInterpolator
lat = lat[::-1]
raw_u = raw_u[:, ::-1, :]
raw_v = raw_v[:, ::-1, :]

# --- load LLC270 grid ---
g  = nc.Dataset(GRID_PATH)
XC = np.array(g.variables['XC'][:])  # (270, 270) lon, -180 to 180
YC = np.array(g.variables['YC'][:])  # (270, 270) lat, ~67.5 to 90
CS = np.array(g.variables['CS'][:])  # (270, 270) cos of rotation angle
SN = np.array(g.variables['SN'][:])  # (270, 270) sin of rotation angle
g.close()
print(f'Grid loaded: XC range=[{XC.min():.1f}, {XC.max():.1f}]  YC range=[{YC.min():.1f}, {YC.max():.1f}]')

# ERA5 lon is 0-360, XC is -180 to 180 -- convert for query
XC_query = np.where(XC < 0, XC + 360.0, XC)
pts = np.column_stack([YC.ravel(), XC_query.ravel()])  # (270*270, 2) as (lat, lon)

# --- build 12-month climatology and interpolate ---
taux_out = np.zeros((12, 270, 270), dtype=np.float64)
tauy_out = np.zeros((12, 270, 270), dtype=np.float64)

for m in range(12):
    idx    = np.arange(m, 360, 12)          # 30 years of this month
    u_clim = raw_u[idx].mean(axis=0)        # (721, 1440) geographic eastward
    v_clim = raw_v[idx].mean(axis=0)        # (721, 1440) geographic northward

    fu = RegularGridInterpolator((lat, lon), u_clim,
                                  method='linear', bounds_error=False, fill_value=None)
    fv = RegularGridInterpolator((lat, lon), v_clim,
                                  method='linear', bounds_error=False, fill_value=None)

    tx_geo = fu(pts).reshape(270, 270)      # geographic eastward on LLC grid
    ty_geo = fv(pts).reshape(270, 270)      # geographic northward on LLC grid

    # rotate geographic -> grid-relative
    # CS = cos(theta), SN = sin(theta), theta = angle from grid-x to east
    tx_grid =  tx_geo * CS + ty_geo * SN
    ty_grid = -tx_geo * SN + ty_geo * CS

    taux_out[m] = tx_grid
    tauy_out[m] = ty_grid

    print(f'  month {m+1:02d}: geo tx=[{tx_geo.min():.4f},{tx_geo.max():.4f}]  '
          f'grid tx=[{tx_grid.min():.4f},{tx_grid.max():.4f}]')

# --- write big-endian float64 binaries ---
taux_path = f'{OUT_DIR}/taux_monthly_rot.bin'
tauy_path = f'{OUT_DIR}/tauy_monthly_rot.bin'

with open(taux_path, 'wb') as f:
    for m in range(12):
        taux_out[m].astype('>f8').tofile(f)

with open(tauy_path, 'wb') as f:
    for m in range(12):
        tauy_out[m].astype('>f8').tofile(f)

expected = 12 * 270 * 270 * 8
for path in [taux_path, tauy_path]:
    sz = os.path.getsize(path)
    ok = 'OK' if sz == expected else f'MISMATCH expected {expected}'
    print(f'{path}: {sz} bytes [{ok}]')

# --- verification plot ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

proj     = ccrs.NorthPolarStereo()
data_crs = ccrs.PlateCarree()
months   = ['Jan','Feb','Mar','Apr','May','Jun',
            'Jul','Aug','Sep','Oct','Nov','Dec']
S = 12

fig, axes = plt.subplots(3, 4, figsize=(18, 13),
                          subplot_kw={'projection': proj})

for m, ax in enumerate(axes.flat):
    ax.set_extent([-180, 180, 62, 90], crs=data_crs)
    ax.add_feature(cfeature.LAND,      facecolor='lightgray', zorder=4)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3,         zorder=5)
    ax.gridlines(linewidth=0.2, color='gray', alpha=0.4)

    speed = np.sqrt(taux_out[m]**2 + tauy_out[m]**2)
    ax.pcolormesh(XC, YC, speed, transform=data_crs,
                  cmap='YlOrRd', vmin=0, vmax=0.15, zorder=2)

    # back-rotate grid-relative -> geographic east/north for display
    tx_disp =  taux_out[m] * CS - tauy_out[m] * SN
    ty_disp =  taux_out[m] * SN + tauy_out[m] * CS

    ax.quiver(XC[::S, ::S], YC[::S, ::S],
              tx_disp[::S, ::S], ty_disp[::S, ::S],
              transform=data_crs,
              scale=1.2, width=0.003, color='black', zorder=3)
    ax.set_title(months[m], fontsize=9)

plt.suptitle('Wind stress -- grid-relative magnitude, geographic vectors', fontsize=12)
plt.tight_layout()
out_fig = f'{OUT_DIR}/wind_verification.png'
plt.savefig(out_fig, dpi=130, bbox_inches='tight')
print(f'Saved verification plot: {out_fig}')
