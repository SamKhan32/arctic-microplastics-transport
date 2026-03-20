import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ERA5_PATH = "original_data/" + 'era5_stress_monthly_1991_2020.nc'
GRID_PATH = "original_data/" + 'ECCO-GRID_06.nc'
OUT_DIR   = "inputs/"

os.makedirs(OUT_DIR, exist_ok=True)

# Load ERA5 -- shape (360, 721, 1440)
ds_era5  = nc.Dataset(ERA5_PATH)
lat_era5 = np.array(ds_era5.variables['latitude'][:])
lon_era5 = np.array(ds_era5.variables['longitude'][:])
iews     = np.array(ds_era5.variables['avg_iews'][:])
inss     = np.array(ds_era5.variables['avg_inss'][:])
ds_era5.close()

if lat_era5[0] > lat_era5[-1]:
    lat_era5 = lat_era5[::-1]
    iews     = iews[:, ::-1, :]
    inss     = inss[:, ::-1, :]

ds_grid = nc.Dataset(GRID_PATH)
xc = np.array(ds_grid.variables['XC'][:])
yc = np.array(ds_grid.variables['YC'][:])
CS = np.array(ds_grid.variables['CS'][:].squeeze())
SN = np.array(ds_grid.variables['SN'][:].squeeze())
ds_grid.close()

xc_query = np.where(xc < 0, xc + 360, xc)
points   = np.column_stack([yc.ravel(), xc_query.ravel()])

taux_stack    = np.zeros((12, 270, 270), dtype=np.float64)
tauy_stack    = np.zeros((12, 270, 270), dtype=np.float64)
taux_geo_save = np.zeros((12, 270, 270), dtype=np.float64)
tauy_geo_save = np.zeros((12, 270, 270), dtype=np.float64)

for m in range(12):
    month_indices = np.arange(m, 360, 12)
    taux_clim     = iews[month_indices].mean(axis=0)
    tauy_clim     = inss[month_indices].mean(axis=0)

    interp_taux = RegularGridInterpolator(
        (lat_era5, lon_era5), taux_clim,
        method='linear', bounds_error=False, fill_value=None)
    interp_tauy = RegularGridInterpolator(
        (lat_era5, lon_era5), tauy_clim,
        method='linear', bounds_error=False, fill_value=None)

    tx_geo = interp_taux(points).reshape(270, 270)
    ty_geo = interp_tauy(points).reshape(270, 270)

    tx_grid =  tx_geo * CS + ty_geo * SN
    ty_grid = -tx_geo * SN + ty_geo * CS

    taux_stack[m]    = tx_grid
    tauy_stack[m]    = ty_grid
    taux_geo_save[m] = tx_geo
    tauy_geo_save[m] = ty_geo

    print(f"Month {m+1:02d}: taux_rot [{tx_grid.min():.4f}, {tx_grid.max():.4f}]  "
          f"tauy_rot [{ty_grid.min():.4f}, {ty_grid.max():.4f}]")

with open(f'{OUT_DIR}/taux_monthly_rot.bin', 'wb') as f:
    for m in range(12):
        taux_stack[m].astype('>f8').tofile(f)

with open(f'{OUT_DIR}/tauy_monthly_rot.bin', 'wb') as f:
    for m in range(12):
        tauy_stack[m].astype('>f8').tofile(f)

expected = 12 * 270 * 270 * 8
for name in ['taux_monthly_rot.bin', 'tauy_monthly_rot.bin']:
    size = os.path.getsize(f'{OUT_DIR}/{name}')
    status = 'OK' if size == expected else f'MISMATCH -- expected {expected}'
    print(f"{name}: {size} bytes [{status}]")

# --- plots ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

month_names = ['Jan','Feb','Mar','Apr','May','Jun',
               'Jul','Aug','Sep','Oct','Nov','Dec']
plot_months = [0, 3, 6, 9]  # Jan Apr Jul Oct
proj     = ccrs.NorthPolarStereo()
data_crs = ccrs.PlateCarree()
s = 12

for label, tx_all, ty_all, fname in [
        ('Pre-rotation (geographic)',  taux_geo_save, tauy_geo_save, 'wind_geo.png'),
        ('Post-rotation (grid-relative)', taux_stack, tauy_stack,   'wind_rot.png')]:

    fig, axes = plt.subplots(2, 2, figsize=(14, 12),
                             subplot_kw={'projection': proj})
    fig.suptitle(f'ERA5 wind stress -- {label}', fontsize=13)

    for ax, m in zip(axes.flat, plot_months):
        tx = tx_all[m]
        ty = ty_all[m]

        # for display, back-rotate grid-relative to geographic east/north
        if 'grid' in label:
            tx_disp =  tx * CS - ty * SN
            ty_disp =  tx * SN + ty * CS
        else:
            tx_disp = tx
            ty_disp = ty

        mag = np.sqrt(tx**2 + ty**2)
        ax.pcolormesh(xc, yc, mag, transform=data_crs,
                      cmap='viridis', vmin=0, vmax=0.12)
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=3)
        ax.coastlines(zorder=4)
        ax.set_extent([-180, 180, 65, 90], crs=data_crs)
        ax.quiver(xc[::s,::s], yc[::s,::s],
                  tx_disp[::s,::s], ty_disp[::s,::s],
                  transform=data_crs, scale=2.0, width=0.004, color='white')
        ax.set_title(month_names[m])
        ax.gridlines(linewidth=0.4, color='gray', alpha=0.5)

    plt.tight_layout()
    plt.savefig(fname, dpi=130)
    print(f'Saved {fname}')