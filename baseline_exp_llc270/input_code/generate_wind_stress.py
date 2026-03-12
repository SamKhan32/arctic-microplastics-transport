import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os

ERA5_PATH  = 'original_data/era5_stress_monthly_1991_2020.nc'
GRID_PATH  = 'original_data/ECCO-GRID_06.nc'
OUT_DIR    = 'mitgcm_input'

os.makedirs(OUT_DIR, exist_ok=True)

# Load ERA5
ds_era5  = nc.Dataset(ERA5_PATH)
lat_era5 = np.array(ds_era5.variables['latitude'][:])
lon_era5 = np.array(ds_era5.variables['longitude'][:])
iews     = np.array(ds_era5.variables['avg_iews'][:])  # (360, 721, 1440)
inss     = np.array(ds_era5.variables['avg_inss'][:])  # (360, 721, 1440)
ds_era5.close()

# January climatology (months 0, 12, 24, ...)
jan_indices = np.arange(0, 360, 12)
taux_clim   = iews[jan_indices].mean(axis=0)  # (721, 1440)
tauy_clim   = inss[jan_indices].mean(axis=0)

# ERA5 latitude is descending (90 to -90); RegularGridInterpolator needs ascending
if lat_era5[0] > lat_era5[-1]:
    lat_era5  = lat_era5[::-1]
    taux_clim = taux_clim[::-1, :]
    tauy_clim = tauy_clim[::-1, :]

interp_taux = RegularGridInterpolator(
    (lat_era5, lon_era5), taux_clim, method='linear', bounds_error=False, fill_value=None)
interp_tauy = RegularGridInterpolator(
    (lat_era5, lon_era5), tauy_clim, method='linear', bounds_error=False, fill_value=None)

# Load LLC270 grid
ds_grid = nc.Dataset(GRID_PATH)
xc = np.array(ds_grid.variables['XC'][:])  # (270, 270), -180 to 180
yc = np.array(ds_grid.variables['YC'][:])  # (270, 270)
ds_grid.close()

# ERA5 lon is 0 to 360; XC is -180 to 180
xc_query = np.where(xc < 0, xc + 360, xc)
points   = np.column_stack([yc.ravel(), xc_query.ravel()])

taux_llc = interp_taux(points).reshape(270, 270)
tauy_llc = interp_tauy(points).reshape(270, 270)

taux_llc.astype('>f8').tofile(f'{OUT_DIR}/taux_jan.bin')
tauy_llc.astype('>f8').tofile(f'{OUT_DIR}/tauy_jan.bin')

print(f"taux range: {taux_llc.min():.4f} to {taux_llc.max():.4f}")
print(f"tauy range: {tauy_llc.min():.4f} to {tauy_llc.max():.4f}")
print("Wrote taux_jan.bin, tauy_jan.bin")
