import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os

ERA5_PATH = 'original_data/era5_stress_monthly_1991_2020.nc'
GRID_PATH = 'original_data/ECCO-GRID_06.nc'
OUT_DIR   = 'mitgcm_input'

os.makedirs(OUT_DIR, exist_ok=True)

# Load ERA5 -- shape (360, 721, 1440)
ds_era5  = nc.Dataset(ERA5_PATH)
lat_era5 = np.array(ds_era5.variables['latitude'][:])
lon_era5 = np.array(ds_era5.variables['longitude'][:])
iews     = np.array(ds_era5.variables['avg_iews'][:])  # (360, 721, 1440)
inss     = np.array(ds_era5.variables['avg_inss'][:])  # (360, 721, 1440)
ds_era5.close()

# ERA5 latitude is descending; flip to ascending for RegularGridInterpolator
if lat_era5[0] > lat_era5[-1]:
    lat_era5 = lat_era5[::-1]
    iews     = iews[:, ::-1, :]
    inss     = inss[:, ::-1, :]

# Load LLC270 grid coordinates and rotation angles
ds_grid = nc.Dataset(GRID_PATH)
xc = np.array(ds_grid.variables['XC'][:])          # (270, 270)
yc = np.array(ds_grid.variables['YC'][:])          # (270, 270)
CS = np.array(ds_grid.variables['CS'][:].squeeze()) # (270, 270)
SN = np.array(ds_grid.variables['SN'][:].squeeze()) # (270, 270)
ds_grid.close()

# ERA5 lon is 0 to 360; XC is -180 to 180
xc_query = np.where(xc < 0, xc + 360, xc)
points   = np.column_stack([yc.ravel(), xc_query.ravel()])

taux_stack = np.zeros((12, 270, 270), dtype=np.float64)
tauy_stack = np.zeros((12, 270, 270), dtype=np.float64)

for m in range(12):
    month_indices = np.arange(m, 360, 12)
    taux_clim     = iews[month_indices].mean(axis=0)  # (721, 1440)
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

    taux_stack[m] = tx_grid
    tauy_stack[m] = ty_grid

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