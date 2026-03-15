# shared/preprocessing/identify_obcs_boundary.py
# Identify open boundary candidate cells on west and east edges of LLC270 Arctic cap tile.
# Run from shared/

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- paths (relative to shared/) ---
BATHY_FILE = 'inputs/bathy_arctic.bin'
GRID_NC    = 'original_data/ECCO-GRID_06.nc'
OUT_FIG    = '../exp3_openbounds/analysis/figures/prerun/southern_boundary_candidates.png'

NX, NY = 270, 270

# --- load bathy and grid ---
bathy = np.fromfile(BATHY_FILE, dtype='>f8').reshape(NY, NX)
g     = nc.Dataset(GRID_NC)
XC    = np.array(g.variables['XC']).squeeze()
YC    = np.array(g.variables['YC']).squeeze()

# --- find wet cells on each edge ---
west_wet = np.where(bathy[:, 0]   < 0)[0]
east_wet = np.where(bathy[:, -1]  < 0)[0]
south_wet = np.where(bathy[0, :]  < 0)[0]
north_wet = np.where(bathy[-1, :] < 0)[0]

# --- print table ---
for label, idxs, get_lon, get_lat in [
    ('West  (i=0,   OB_Iwest)',  west_wet,  lambda j: XC[j, 0],    lambda j: YC[j, 0]),
    ('East  (i=269, OB_Ieast)',  east_wet,  lambda j: XC[j, -1],   lambda j: YC[j, -1]),
    ('South (j=0,   OB_Jsouth)', south_wet, lambda i: XC[0, i],    lambda i: YC[0, i]),
    ('North (j=269, OB_Jnorth)', north_wet, lambda i: XC[-1, i],   lambda i: YC[-1, i]),
]:
    print(f"\n--- {label}: {len(idxs)} wet cells ---")
    print(f"  {'idx_py':>8} {'idx_f90':>8} {'depth_m':>10} {'lon':>8} {'lat':>8}")
    for idx in idxs:
        print(f"  {idx:>8} {idx+1:>8} {bathy[idx,0] if 'West' in label or 'East' in label else bathy[0,idx]:>10.1f}"
              f" {get_lon(idx):>8.2f} {get_lat(idx):>8.2f}")

# --- classify west edge cells ---
# Fram Strait: lon 0-20E, lat ~72N  -> j~110-145
# Barents:     lon 30-55E           -> j~15-80
# E Greenland: lon ~-24 to -10      -> j~200-212
fram_west   = west_wet[(XC[west_wet, 0] >= 0)  & (XC[west_wet, 0] <= 20)]
barents     = west_wet[(XC[west_wet, 0] >  20)]
e_greenland = west_wet[(XC[west_wet, 0] <  0)]

# --- classify east edge cells ---
# Bering/Chukchi: lon 160E to 180 and -180 to -160
# Beaufort shelf: lon < -130
bering      = east_wet[(XC[east_wet, -1] >= 150) | (XC[east_wet, -1] <= -160)]
beaufort    = east_wet[(XC[east_wet, -1] > -160) & (XC[east_wet, -1] < -130)]

print(f"\n--- Classification summary ---")
print(f"  West/Fram Strait (lon 0-20E):    {len(fram_west)} cells, j_f90 {fram_west[0]+1}-{fram_west[-1]+1}")
print(f"  West/Barents (lon >20E):         {len(barents)} cells, j_f90 {barents[0]+1}-{barents[-1]+1}")
print(f"  West/E.Greenland (lon <0):       {len(e_greenland)} cells, j_f90 {e_greenland[0]+1}-{e_greenland[-1]+1}")
print(f"  East/Bering-Chukchi:             {len(bering)} cells, j_f90 {bering[0]+1}-{bering[-1]+1}")
print(f"  East/Beaufort shelf:             {len(beaufort)} cells, j_f90 {beaufort[0]+1}-{beaufort[-1]+1}")

# --- map ---
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.NorthPolarStereo()}, figsize=(9, 9))
ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.gridlines(draw_labels=False)

# all wet cells faintly
wet_all = np.where(bathy < 0)
ax.scatter(XC[wet_all], YC[wet_all], s=0.3, c='steelblue',
           transform=ccrs.PlateCarree(), alpha=0.15)

# boundary cells colored by region
def plot_edge(idxs, get_lon, get_lat, color, label):
    lons = [get_lon(i) for i in idxs]
    lats = [get_lat(i) for i in idxs]
    ax.scatter(lons, lats, s=12, c=color, transform=ccrs.PlateCarree(),
               label=label, zorder=5)

plot_edge(fram_west,   lambda j: XC[j,0],  lambda j: YC[j,0],  'red',    'Fram Strait (W edge)')
plot_edge(barents,     lambda j: XC[j,0],  lambda j: YC[j,0],  'orange', 'Barents (W edge)')
plot_edge(e_greenland, lambda j: XC[j,0],  lambda j: YC[j,0],  'purple', 'E Greenland (W edge)')
plot_edge(bering,      lambda j: XC[j,-1], lambda j: YC[j,-1], 'green',  'Bering/Chukchi (E edge)')
plot_edge(beaufort,    lambda j: XC[j,-1], lambda j: YC[j,-1], 'cyan',   'Beaufort shelf (E edge)')

ax.legend(loc='lower left', fontsize=8)
ax.set_title('LLC270 Arctic cap -- OBC candidate cells by region')
plt.tight_layout()
plt.savefig(OUT_FIG, dpi=150)
print(f"\nFigure saved to {OUT_FIG}")