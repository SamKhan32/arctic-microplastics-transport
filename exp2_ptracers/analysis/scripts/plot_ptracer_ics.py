# plot_ptracer_ics.py
# Plots both PTRACER IC fields on polar stereographic projection for placement verification
# Run from exp2_ptracers/

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

GRID_NC = "../shared/original_data/ECCO-GRID_06.nc"
TR1     = "run/input/ptracers_ic_tr1.bin"
TR2     = "run/input/ptracers_ic_tr2.bin"

with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)

tr1 = np.fromfile(TR1, dtype='>f8').reshape(1, 270, 270)[0]
tr2 = np.fromfile(TR2, dtype='>f8').reshape(1, 270, 270)[0]

fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                         subplot_kw={'projection': ccrs.NorthPolarStereo()})

titles  = ["PTRACER01 -- Real samples (n=68)", "PTRACER02 -- Gyre/drift placeholders"]
tracers = [tr1, tr2]

for ax, tracer, title in zip(axes, tracers, titles):
    ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.gridlines(draw_labels=False, linewidth=0.4, color='gray', linestyle='--')

    masked = np.ma.masked_where(tracer == 0, tracer)
    ax.pcolormesh(XC, YC, masked, transform=ccrs.PlateCarree(),
                  cmap='plasma', shading='auto')

    ax.set_title(title)

plt.tight_layout()
plt.savefig("analysis/figures/prerun/ptracer_ic_placement.png", dpi=150)
print("Saved to analysis/figures/prerun/ptracer_ic_placement.png")