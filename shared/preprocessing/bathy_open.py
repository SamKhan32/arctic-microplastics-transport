"""
bathy_open.py
-------------
Interpolates IBCAO v5 400m bathymetry onto the LLC270 Arctic cap
curvilinear grid and writes a MITgcm-ready binary for open-boundary
experiments. Unlike the closed-basin version, Bering Strait is left
open and no connected-region filtering is applied -- edge cells retain
their true depths so OBCS boundary cells are wet.

Run from shared/
"""

import numpy as np
import rasterio
from pyproj import Transformer
from scipy.ndimage import map_coordinates
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- paths (relative to shared/) ---
TIF_PATH   = 'original_data/ibcao_v5_1_2025_400m_ice.tif'
GRID_NC    = 'original_data/ECCO-GRID_06.nc'
OUT_BIN    = 'inputs/bathy_arctic_open.bin'
OUT_PLOT   = '../exp3_openbounds/analysis/figures/prerun/bathy_open_check.png'

NX, NY = 270, 270


def load_grid(grid_nc):
    g   = nc.Dataset(grid_nc)
    XC  = np.array(g.variables['XC']).squeeze()  # (NY, NX)
    YC  = np.array(g.variables['YC']).squeeze()
    print(f"Grid loaded: XC shape={XC.shape}, lon [{XC.min():.1f}, {XC.max():.1f}]")
    return XC, YC


def interpolate_bathy(tif_path, LON, LAT):
    with rasterio.open(tif_path) as src:
        print(f"IBCAO CRS:        {src.crs}")
        print(f"IBCAO resolution: {src.res} m")
        print(f"IBCAO shape:      {src.shape}")
        data      = src.read(1).astype('float32')
        transform = src.transform
        ibcao_crs = src.crs
        if src.nodata is not None:
            data = np.where(data == src.nodata, 0.0, data)

    # reproject grid centers to IBCAO polar stereographic
    transformer = Transformer.from_crs(
        'EPSG:4326',
        ibcao_crs.to_epsg() if ibcao_crs.to_epsg() else ibcao_crs.to_wkt(),
        always_xy=True
    )
    x_ps, y_ps = transformer.transform(LON, LAT)

    col_float = (x_ps - transform.c) / transform.a
    row_float = (y_ps - transform.f) / transform.e

    print(f"Pixel col range: {col_float.min():.1f} to {col_float.max():.1f}")
    print(f"Pixel row range: {row_float.min():.1f} to {row_float.max():.1f}")

    print("Interpolating (bilinear)...")
    bathy = map_coordinates(
        data,
        [row_float, col_float],
        order=1,
        mode='nearest'
    )

    # mask land: depth >= 0 -> 0.0
    # no Bering closure, no connected-region filtering
    bathy = np.where(bathy >= 0, 0.0, bathy)

    print(f"Bathy stats:")
    print(f"  min: {bathy.min():.1f} m  max: {bathy.max():.1f} m")
    print(f"  ocean cells: {(bathy < 0).sum()}  land cells: {(bathy == 0).sum()}")
    print(f"  West edge wet (i=0):   {(bathy[:, 0] < 0).sum()}")
    print(f"  East edge wet (i=269): {(bathy[:, -1] < 0).sum()}")

    return bathy.astype('>f8')


def write_binary(bathy, out_path):
    bathy.flatten(order='C').tofile(out_path)
    print(f"Written: {out_path}  dtype={bathy.dtype}  shape={bathy.shape}")


def plot_bathy(bathy, LON, LAT, out_plot):
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.NorthPolarStereo()}, figsize=(9, 9))
    ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
    depth = np.where(bathy < 0, bathy, np.nan)
    ax.pcolormesh(LON, LAT, depth, transform=ccrs.PlateCarree(),
                  cmap='Blues_r', shading='auto')
    ax.add_feature(cfeature.LAND, color='lightgray', zorder=4)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=5)
    ax.gridlines(linewidth=0.4, color='gray', alpha=0.5, linestyle='--')
    ax.set_title('LLC270 Arctic cap -- open boundary bathy\n(Bering open, no region filtering)')
    plt.tight_layout()
    plt.savefig(out_plot, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {out_plot}")


if __name__ == '__main__':
    LON, LAT = load_grid(GRID_NC)
    bathy    = interpolate_bathy(TIF_PATH, LON, LAT)
    write_binary(bathy, OUT_BIN)
    plot_bathy(bathy, LON, LAT, OUT_PLOT)