import numpy as np
import netCDF4 as nc

g = nc.Dataset('original_data/ECCO-GRID_06.nc')
CS = np.array(g['CS'][:].squeeze())
SN = np.array(g['SN'][:].squeeze())
g.close()

tx_geo = np.fromfile('mitgcm_input/taux_jan.bin', dtype='>f4').reshape(270, 270)
ty_geo = np.fromfile('mitgcm_input/tauy_jan.bin', dtype='>f4').reshape(270, 270)

tx_grid =  tx_geo * CS + ty_geo * SN
ty_grid = -tx_geo * SN + ty_geo * CS

tx_grid.astype('>f4').tofile('mitgcm_input/taux_jan_rot.bin')
ty_grid.astype('>f4').tofile('mitgcm_input/tauy_jan_rot.bin')

print('taux_rot: min=%g max=%g' % (tx_grid.min(), tx_grid.max()))
print('tauy_rot: min=%g max=%g' % (ty_grid.min(), ty_grid.max()))

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

g = nc.Dataset('original_data/ECCO-GRID_06.nc')
XC = np.array(g['XC'][:].squeeze())
YC = np.array(g['YC'][:].squeeze())
g.close()

tx_geo = np.fromfile('mitgcm_input/taux_jan.bin', dtype='>f4').reshape(270, 270)
ty_geo = np.fromfile('mitgcm_input/tauy_jan.bin', dtype='>f4').reshape(270, 270)
tx_rot = np.fromfile('mitgcm_input/taux_jan_rot.bin', dtype='>f4').reshape(270, 270)
ty_rot = np.fromfile('mitgcm_input/tauy_jan_rot.bin', dtype='>f4').reshape(270, 270)

proj = ccrs.NorthPolarStereo()
data_crs = ccrs.PlateCarree()

step = 10
fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                          subplot_kw={'projection': proj})

for ax, tx, ty, title in [
    (axes[0], tx_geo, ty_geo, 'Geographic (original)'),
    (axes[1], tx_rot, ty_rot, 'Grid-relative (rotated)')
]:
    ax.set_extent([-180, 180, 60, 90], crs=data_crs)
    ax.add_feature(cfeature.LAND, zorder=1)
    ax.add_feature(cfeature.COASTLINE, zorder=2)
    ax.gridlines()
    mag = np.sqrt(tx**2 + ty**2)
    im = ax.pcolormesh(XC, YC, mag, transform=data_crs,
                       cmap='viridis', vmin=0, vmax=0.3)
    ax.quiver(XC[::step, ::step], YC[::step, ::step],
              tx[::step, ::step], ty[::step, ::step],
              transform=data_crs, scale=4, color='white', width=0.003)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, shrink=0.7, label='N/m²')

plt.suptitle('ERA5 January Wind Stress — Arctic Cap')
plt.tight_layout()
plt.savefig('wind_stress_check.png', dpi=150)
plt.show()