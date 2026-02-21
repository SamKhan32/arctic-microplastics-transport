"""
bathy.py
--------
Interpolates IBCAO v5 400m bathymetry onto the BASELINE_EXP
curvilinear grid and writes a MITgcm-ready binary file.

Workflow:
  1. Load curvilinear grid cell centers (LON, LAT) from config
  2. Reproject grid centers into IBCAO polar stereographic coords
  3. Sample IBCAO data at those coords via bilinear interpolation
  4. Mask land (depth >= 0) → 0.0 (MITgcm convention)
  5. Write big-endian float32 binary in Fortran order
  6. Plot to verify
"""

import numpy as np
import rasterio
from pyproj import Transformer
from scipy.ndimage import map_coordinates
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
from pathlib import Path

# Import shared config — config.py sits one level up in BASELINE_EXP/
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import (
    TIF_PATH, MITGCM_INPUT_DIR,
    NX, NY, make_grid
)


def interpolate_bathy(tif_path, LON, LAT, dtype='>f4'):
    """
    Interpolate IBCAO GeoTIFF onto curvilinear grid cell centers.

    Parameters
    ----------
    tif_path : Path
        Path to IBCAO GeoTIFF
    LON, LAT : np.ndarray, shape (NY, NX)
        Geographic coordinates of curvilinear grid cell centers
    dtype : str
        MITgcm binary dtype, default big-endian float32

    Returns
    -------
    bathy : np.ndarray, shape (NY, NX)
        Bathymetry in meters, land = 0.0, ocean = negative depth
    """
    with rasterio.open(tif_path) as src:
        print(f"IBCAO CRS:        {src.crs}")
        print(f"IBCAO resolution: {src.res} m")
        print(f"IBCAO shape:      {src.shape}")
        print(f"IBCAO bounds:     {src.bounds}")
        print(f"IBCAO nodata:     {src.nodata}")

        data = src.read(1).astype('float32')  # (nrows, ncols)
        transform = src.transform
        ibcao_crs = src.crs

        # Replace nodata with land (0.0)
        if src.nodata is not None:
            data = np.where(data == src.nodata, 0.0, data)

    # Step 1: reproject curvilinear grid centers (lat/lon) →
    #         IBCAO polar stereographic (meters)
    print("\nReprojecting grid centers to IBCAO CRS...")
    transformer = Transformer.from_crs(
        "EPSG:4326",       # geographic lat/lon
        ibcao_crs.to_epsg() if ibcao_crs.to_epsg() else ibcao_crs.to_wkt(),
        always_xy=True
    )
    x_ps, y_ps = transformer.transform(LON, LAT)

    # Step 2: convert projected coords → pixel (row, col) in IBCAO array
    # rasterio transform: x = T.c + col*T.a,  y = T.f + row*T.e
    col_float = (x_ps - transform.c) / transform.a
    row_float = (y_ps - transform.f) / transform.e

    print(f"Pixel col range: {col_float.min():.1f} to {col_float.max():.1f}")
    print(f"Pixel row range: {row_float.min():.1f} to {row_float.max():.1f}")
    print(f"IBCAO array size: {data.shape[0]} rows x {data.shape[1]} cols")

    # Step 3: bilinear interpolation (order=1)
    # map_coordinates expects [row_coords, col_coords]
    # mode='nearest' handles edge cells that fall just outside bounds
    print("Interpolating (bilinear)...")
    bathy = map_coordinates(
        data,
        [row_float, col_float],
        order=1,
        mode='nearest'
    )
    
    # Step 4: mask land — anything >= 0 is land or sea ice surface → 0.0
    bathy = np.where(bathy >= 0, 0.0, bathy)
    from scipy.ndimage import label
    # Close Bering Strait
    # Geographic location: roughly 65-66°N, 168-170°W
    bering_mask = (
        (LAT > 64.5) & (LAT < 66.5) &
        (LON > -171.0) & (LON < -166.0) &
        (bathy < 0)
    )
    print(f"Bering Strait cells to close: {bering_mask.sum()}")
    bathy = np.where(bering_mask, 0.0, bathy)
    # Find all connected ocean regions
    ocean_mask = bathy < 0
    labeled, num_features = label(ocean_mask)
    print(f"Found {num_features} disconnected ocean regions")

    # Find the largest connected region (the Arctic Ocean)
    region_sizes = np.bincount(labeled.ravel())
    region_sizes[0] = 0  # ignore land (label 0)
    largest_region = region_sizes.argmax()

    # Keep only the largest region, set everything else to land
    bathy = np.where(labeled == largest_region, bathy, 0.0)

    print(f"Removed {num_features - 1} isolated ocean regions")
    print(f"\nInterpolated bathy stats:")
    print(f"  min depth: {bathy.min():.1f} m")
    print(f"  max depth: {bathy.max():.1f} m")
    print(f"  land cells: {(bathy == 0).sum()} / {bathy.size}")
    print(f"  ocean cells: {(bathy < 0).sum()} / {bathy.size}")

    return bathy.astype(dtype)


def write_mitgcm_binary(bathy, out_path, dtype='>f4'):
    """
    Write bathymetry as MITgcm-ready binary.
    MITgcm reads (NX, NY) with Fortran column-major order.
    bathy input shape: (NY, NX)
    """
    bathy_T = bathy.T  # → (NX, NY)
    bathy_T.flatten(order='F').tofile(out_path)
    print(f"\nWritten: {out_path}")
    print(f"  Shape written (NX, NY): {bathy_T.shape}")
    print(f"  dtype: {dtype}")


def plot_bathy(bathy, LON, LAT):
    """Quick verification plot on Arctic polar stereographic projection."""
    fig = plt.figure(figsize=(10, 10))
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

    # Plot bathymetry
    # bathy shape (NY, NX), LON/LAT shape (NY, NX)
    depth = np.where(bathy < 0, bathy, np.nan)  # mask land for colormap
    import matplotlib.colors as mcolors

    # Define depth boundaries and colors
    bounds = [0, -50, -200, -500, -1000, -2000, -3500, -5500]
    colors = [
        '#d4eaf7',  # 0 to -50m      — very shallow shelf
        '#a8d5f0',  # -50 to -200m   — shallow shelf
        '#6ab4e8',  # -200 to -500m  — shelf break
        '#3d8ec4',  # -500 to -1000m — upper slope
        '#1f5f8b',  # -1000 to -2000m
        '#0d3a5c',  # -2000 to -3500m
        '#051d30',  # -3500 to -5500m — deep basin
    ]

    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds[::-1], cmap.N)  # reverse since depths are negative

    im = ax.pcolormesh(
        LON, LAT, depth,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        shading='auto'
    )
    plt.colorbar(im, ax=ax, label='Depth (m)', shrink=0.7)
    # Add land AFTER the pcolormesh so it paints over the edges
    ax.add_feature(cfeature.LAND, color='#c8b89a', zorder=4)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=5)

    ax.add_feature(cfeature.LAND, color='#c8b89a', zorder=3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=4)
    ax.gridlines(linewidth=0.4, color='gray', alpha=0.5, linestyle='--')

    ax.set_title(
        f'BASELINE_EXP Bathymetry\n'
        f'IBCAO v5 interpolated onto {NX}×{NY} curvilinear grid',
        fontsize=12
    )
    plt.tight_layout()
    plt.savefig(
        MITGCM_INPUT_DIR.parent / 'analysis' / 'bathy_check.png',
        dpi=150, bbox_inches='tight'
    )
    print("Plot saved to analysis/bathy_check.png")
    plt.show()


if __name__ == "__main__":
    # Build curvilinear grid
    print("Building curvilinear grid...")
    LON, LAT = make_grid()
    print(f"Grid shape: {LON.shape} (NY={NY}, NX={NX})")

    # Interpolate IBCAO onto grid
    bathy = interpolate_bathy(TIF_PATH, LON, LAT)

    # Write MITgcm binary
    out_file = MITGCM_INPUT_DIR / "bathy.bin"
    write_mitgcm_binary(bathy, out_file)

    # Verify with plot
    plot_bathy(bathy, LON, LAT)