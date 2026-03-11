import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Load ERA5
ds_era5 = nc.Dataset('original_data/era5_stress_monthly_1991_2020.nc')
lat_era5 = np.array(ds_era5.variables['latitude'][:])
lon_era5 = np.array(ds_era5.variables['longitude'][:])
iews = np.array(ds_era5.variables['avg_iews'][:])  # (360, 721, 1440)
inss = np.array(ds_era5.variables['avg_inss'][:])  # (360, 721, 1440)

# January climatology - months 0, 12, 24, ... (0-indexed)
jan_indices = np.arange(0, 360, 12)
taux_clim = iews[jan_indices].mean(axis=0)  # (721, 1440)
tauy_clim = inss[jan_indices].mean(axis=0)  # (721, 1440)

print(f"ERA5 lat range: {lat_era5.min():.1f} to {lat_era5.max():.1f}")
print(f"ERA5 lon range: {lon_era5.min():.1f} to {lon_era5.max():.1f}")
print(f"taux range: {taux_clim.min():.4f} to {taux_clim.max():.4f}")
print(f"tauy range: {tauy_clim.min():.4f} to {tauy_clim.max():.4f}")

# Load LLC270 grid coordinates
ds_grid = nc.Dataset('original_data/ECCO-GRID_06.nc')
xc = np.array(ds_grid.variables['XC'][:])  # (270, 270), -180 to 180
yc = np.array(ds_grid.variables['YC'][:])  # (270, 270), 67.5 to 89.9

# ERA5 latitude is descending (90 to -90), need ascending for interpolator
if lat_era5[0] > lat_era5[-1]:
    lat_era5 = lat_era5[::-1]
    taux_clim = taux_clim[::-1, :]
    tauy_clim = tauy_clim[::-1, :]

# Build interpolators
interp_taux = RegularGridInterpolator(
    (lat_era5, lon_era5), taux_clim, method='linear', bounds_error=False, fill_value=None)
interp_tauy = RegularGridInterpolator(
    (lat_era5, lon_era5), tauy_clim, method='linear', bounds_error=False, fill_value=None)

# Interpolate to LLC270 grid points
# XC is -180 to 180, ERA5 lon is 0 to 360 - need to handle wrap
xc_query = np.where(xc < 0, xc + 360, xc)
points = np.column_stack([yc.ravel(), xc_query.ravel()])

taux_llc = interp_taux(points).reshape(270, 270)
tauy_llc = interp_tauy(points).reshape(270, 270)

print(f"Interpolated taux range: {taux_llc.min():.4f} to {taux_llc.max():.4f}")
print(f"Interpolated tauy range: {tauy_llc.min():.4f} to {tauy_llc.max():.4f}")

# Write binary files
taux_llc.astype('>f4').tofile('original_data/llc270_grid/taux_jan.bin')
tauy_llc.astype('>f4').tofile('original_data/llc270_grid/tauy_jan.bin')

print("Wind stress files written.")

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                         subplot_kw={'projection': ccrs.NorthPolarStereo()})

titles = ['Taux (Eastward Wind Stress)', 'Tauy (Northward Wind Stress)']
fields = [taux_llc, tauy_llc]

for ax, field, title in zip(axes, fields, titles):
    ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())
    p = ax.pcolormesh(xc, yc, field,
                      transform=ccrs.PlateCarree(),
                      cmap='RdBu_r', vmin=-0.5, vmax=0.5)
    ax.coastlines(resolution='50m')
    ax.gridlines()
    plt.colorbar(p, ax=ax, label='N/m²', shrink=0.7)
    ax.set_title(title)

plt.suptitle('ERA5 January Climatology Wind Stress - LLC270 Arctic Cap')
plt.tight_layout()
plt.savefig('wind_stress_check.png', dpi=100)
plt.show()