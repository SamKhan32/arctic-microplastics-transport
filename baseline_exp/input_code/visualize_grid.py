"""
visualize_grid.py
-----------------
Visualizes a displaced-pole curvilinear grid over the Arctic.
Lets you experiment with:
  - Rotated pole location (pole_lon, pole_lat)
  - Grid resolution (dx_km)
  - Grid dimensions (nx, ny)

Requirements:
    pip install numpy matplotlib cartopy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ─────────────────────────────────────────────
# PARAMETERS — tweak these to explore
# ─────────────────────────────────────────────

# Displaced pole location — over Greenland
POLE_LON = -40.0   # degrees East (negative = West)
POLE_LAT =  75.0   # degrees North

# Grid resolution in km (approximate)
DX_KM = 15.0

# Grid dimensions (nx, ny) — number of cells
# At 5km res, ~2000x2000 covers roughly 10,000km x 10,000km
# Start small for visualization, scale up later
NX = 1530
NY = 520

            # accounting for cos factor at mid-rotated-lat)

# Rotated coordinate extent — in the *rotated* system
# Think of these as degrees in the rotated frame
# -90 to 90 in rotated lat covers most of the Arctic
ROT_LON_MIN = -90.0
ROT_LON_MAX =  90
ROT_LAT_MIN =  40.0   # controls how far south — tune this
ROT_LAT_MAX =  90.0   # rotated north pole, leave this

# How many grid lines to show (stride for clarity)
STRIDE = 40

# ─────────────────────────────────────────────
# CORE MATH: rotated pole transformation
# ─────────────────────────────────────────────

def rotate_coords(lon, lat, pole_lon, pole_lat):
    """
    Transform geographic (lon, lat) to rotated-pole coordinates.
    
    The rotated pole is placed at (pole_lon, pole_lat) in geographic coords.
    In the rotated system, this point becomes (0, 90).
    
    Returns rotated (rlon, rlat) in degrees.
    """
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    plon_r = np.deg2rad(pole_lon)
    plat_r = np.deg2rad(pole_lat)

    # Standard rotated-pole forward transform
    sin_plat = np.sin(plat_r)
    cos_plat = np.cos(plat_r)

    sin_lat = np.sin(lat_r)
    cos_lat = np.cos(lat_r)
    cos_dlon = np.cos(lon_r - plon_r)

    # Rotated latitude
    rlat = np.arcsin(
        sin_plat * sin_lat + cos_plat * cos_lat * cos_dlon
    )

    # Rotated longitude
    rlon = np.arctan2(
        cos_lat * np.sin(lon_r - plon_r),
        cos_plat * sin_lat - sin_plat * cos_lat * cos_dlon
    )

    return np.rad2deg(rlon), np.rad2deg(rlat)


def unrotate_coords(rlon, rlat, pole_lon, pole_lat):
    """
    Inverse: rotated-pole (rlon, rlat) → geographic (lon, lat).
    This is what you use to find the true lat/lon of each grid cell.
    """
    rlon_r = np.deg2rad(rlon)
    rlat_r = np.deg2rad(rlat)
    plon_r = np.deg2rad(pole_lon)
    plat_r = np.deg2rad(pole_lat)

    sin_plat = np.sin(plat_r)
    cos_plat = np.cos(plat_r)

    sin_rlat = np.sin(rlat_r)
    cos_rlat = np.cos(rlat_r)
    cos_rlon = np.cos(rlon_r)

    # Geographic latitude
    lat = np.arcsin(
        sin_plat * sin_rlat + cos_plat * cos_rlat * cos_rlon
    )

    # Geographic longitude
    lon = np.arctan2(
        cos_rlat * np.sin(rlon_r),
        cos_plat * sin_rlat - sin_plat * cos_rlat * cos_rlon
    ) + np.deg2rad(pole_lon)

    return np.rad2deg(lon), np.rad2deg(lat)


# ─────────────────────────────────────────────
# BUILD THE GRID
# ─────────────────────────────────────────────

# Define a regular grid in rotated coordinates
rlon_1d = np.linspace(ROT_LON_MIN, ROT_LON_MAX, NX)
rlat_1d = np.linspace(ROT_LAT_MIN, ROT_LAT_MAX, NY)
RLON, RLAT = np.meshgrid(rlon_1d, rlat_1d)  # shape (NY, NX)

# Convert to geographic coordinates — these become your XC, YC arrays
LON, LAT = unrotate_coords(RLON, RLAT, POLE_LON, POLE_LAT)

print(f"Grid shape: {LON.shape}")
print(f"Geographic lon range: {LON.min():.1f} to {LON.max():.1f}")
print(f"Geographic lat range: {LAT.min():.1f} to {LAT.max():.1f}")
print(f"Approximate domain coverage at 5km resolution:")
print(f"  {NX * DX_KM:.0f} km x {NY * DX_KM:.0f} km")

# ─────────────────────────────────────────────
# PLOT
# ─────────────────────────────────────────────
def plot_line_safe(ax, lons, lats, **kwargs):
    """Plot a line, breaking it where it crosses the antimeridian."""
    lons = np.array(lons)
    lats = np.array(lats)
    
    # Find where longitude jumps by more than 180° (antimeridian crossing)
    breaks = np.where(np.abs(np.diff(lons)) > 180)[0] + 1
    
    # Split into segments at those breaks
    segments = np.split(np.arange(len(lons)), breaks)
    
    for seg in segments:
        if len(seg) > 1:
            ax.plot(lons[seg], lats[seg], transform=ccrs.PlateCarree(), **kwargs)
fig = plt.figure(figsize=(12, 10))

# North Polar Stereographic projection — natural view for Arctic
proj = ccrs.NorthPolarStereo(central_longitude=0)
ax = fig.add_subplot(1, 1, 1, projection=proj)

# Set map extent to Arctic
ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

# Background features
ax.add_feature(cfeature.OCEAN, color='#d0e8f5', zorder=0)
ax.add_feature(cfeature.LAND, color='#c8b89a', zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='#555', zorder=2)
ax.add_feature(cfeature.BORDERS, linewidth=0.3, color='#888', zorder=2)

# Draw grid lines (every STRIDE-th row and column)
# i-lines (constant rotated longitude — meridians in rotated system)
for i in range(0, NX, STRIDE):
    plot_line_safe(ax, LON[:, i], LAT[:, i], color='steelblue', linewidth=0.5, alpha=0.7)

for j in range(0, NY, STRIDE):
    plot_line_safe(ax, LON[j, :], LAT[j, :], color='steelblue', linewidth=0.5, alpha=0.7)

# Mark the displaced pole location
ax.plot(POLE_LON, POLE_LAT, 'r*', markersize=15,
        transform=ccrs.PlateCarree(), zorder=5,
        label=f'Displaced pole ({POLE_LON}°E, {POLE_LAT}°N)')

# Geographic grid lines for reference
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
   #               linewidth=0.5, color='gray', alpha=0.4, linestyle='--')

# Arctic Circle reference
arctic_circle_lat = 66.5
theta = np.linspace(0, 2 * np.pi, 360)
ac_lon = np.rad2deg(theta)
ac_lat = np.full_like(ac_lon, arctic_circle_lat)
ax.plot(ac_lon, ac_lat, 'r--', linewidth=1.0, alpha=0.6,
        transform=ccrs.PlateCarree(), label='Arctic Circle (66.5°N)', zorder=4)

ax.legend(loc='lower left', fontsize=9)
ax.set_title(
    f'Displaced-Pole Curvilinear Grid\n'
    f'Pole at ({POLE_LON}°E, {POLE_LAT}°N) | '
    f'{NX}×{NY} cells | ~{DX_KM}km resolution | '
    f'Grid lines every {STRIDE} cells',
    fontsize=11
)

plt.tight_layout()
plt.savefig('displaced_pole_grid.png', dpi=150, bbox_inches='tight')
print("\nSaved to displaced_pole_grid.png")
plt.show()
# Approximate cell size check at a few points
R_earth = 6371000  # meters
dy = np.deg2rad(rlat_1d[1] - rlat_1d[0]) * R_earth
print(f"Approximate dy: {dy/1000:.1f} km")

# Check isotropy at a few geographic locations
# by looking at dx/dy ratio across the grid