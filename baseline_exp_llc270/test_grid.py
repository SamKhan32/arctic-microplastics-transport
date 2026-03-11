import netCDF4 as nc
import numpy as np
import os

path = "original_data/ECCO-GRID_06.nc"
ds = nc.Dataset(path)

os.makedirs('original_data/llc270_grid', exist_ok=True)

def get(varname):
    return np.array(ds.variables[varname][:], dtype=np.float64)

fields = [
    get('XC'),   # 1  XC
    get('YC'),   # 2  YC
    get('dxG'),  # 3  DXF (approximation)
    get('dyG'),  # 4  DYF (approximation)
    get('rA'),   # 5  RAC
    get('XG'),   # 6  XG
    get('YG'),   # 7  YG
    get('dxC'),  # 8  DXV
    get('dyC'),  # 9  DYU
    get('rAz'),  # 10 RAZ
    get('dxC'),  # 11 DXC
    get('dyC'),  # 12 DYC
    get('rAw'),  # 13 RAW
    get('rAs'),  # 14 RAS
    get('dxG'),  # 15 DXG
    get('dyG'),  # 16 DYG
]

with open('original_data/llc270_grid/arctic_cap_v2.mitgrid', 'wb') as f:
    for field in fields:
        field.astype('>f8').tofile(f)

data = np.fromfile('original_data/llc270_grid/arctic_cap_v2.mitgrid', dtype='>f8')
print(f"Elements written: {len(data)}")
print(f"Expected:         {16 * 270 * 270}")
print(f"Match: {len(data) == 16 * 270 * 270}")

depth = np.array(ds.variables['Depth'][:], dtype=np.float64)
# MITgcm convention: ocean negative, land = 0
bathy = np.where(depth > 0, -depth, 0.0)
bathy.astype('>f8').tofile('original_data/llc270_grid/bathy_arctic.bin')
print(f"Bathymetry written: {bathy.shape}, min={bathy.min():.1f}, max={bathy.max():.1f}")

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

xc = np.array(ds.variables['XC'][:])
yc = np.array(ds.variables['YC'][:])
bathy = np.fromfile('original_data/llc270_grid/bathy_arctic.bin', dtype='>f8').reshape(270, 270)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())

p = ax.pcolormesh(xc, yc, bathy,
                  transform=ccrs.PlateCarree(),
                  cmap='Blues_r')

ax.coastlines(resolution='50m')
ax.gridlines()
plt.colorbar(p, label='Depth (m)', shrink=0.7)
plt.title('Arctic Cap Bathymetry (bathy_arctic.bin)')
plt.savefig('arctic_bathy_check.png', dpi=100)
plt.show()