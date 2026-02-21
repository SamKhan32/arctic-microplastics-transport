"""
config.py
---------
Canonical parameters for BASELINE_EXP.
All input_code scripts import from here.
Changing a value here propagates to all scripts automatically.
"""

from pathlib import Path
import numpy as np

# ─────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────
EXP_DIR          = Path(__file__).resolve().parent
TIF_PATH         = EXP_DIR / "original_data" / "ibcao_v5_1_2025_400m_ice.tif"
MITGCM_INPUT_DIR = EXP_DIR / "mitgcm_input"

# ─────────────────────────────────────────────
# CURVILINEAR GRID PARAMETERS
# ─────────────────────────────────────────────
POLE_LON = -40.0   # Displaced pole longitude (over Greenland)
POLE_LAT =  75.0   # Displaced pole latitude  (over Greenland)

NX = 1530          # Grid cells in x (rotated longitude direction)
NY =  520          # Grid cells in y (rotated latitude direction)

ROT_LON_MIN = -90.0
ROT_LON_MAX =  90.0
ROT_LAT_MIN =  40.0
ROT_LAT_MAX =  90.0

# ─────────────────────────────────────────────
# GRID GENERATION (shared function)
# ─────────────────────────────────────────────

def make_grid():
    """
    Build curvilinear grid cell centers (LON, LAT) in geographic coordinates.
    Returns LON, LAT as 2D arrays of shape (NY, NX).
    """
    rlon_1d = np.linspace(ROT_LON_MIN, ROT_LON_MAX, NX)
    rlat_1d = np.linspace(ROT_LAT_MIN, ROT_LAT_MAX, NY)
    RLON, RLAT = np.meshgrid(rlon_1d, rlat_1d)  # (NY, NX)
    LON, LAT = unrotate_coords(RLON, RLAT, POLE_LON, POLE_LAT)
    return LON, LAT


def unrotate_coords(rlon, rlat, pole_lon, pole_lat):
    """
    Rotated-pole (rlon, rlat) → geographic (lon, lat).
    All inputs/outputs in degrees.
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

    lat = np.arcsin(
        sin_plat * sin_rlat + cos_plat * cos_rlat * cos_rlon
    )
    lon = np.arctan2(
        cos_rlat * np.sin(rlon_r),
        cos_plat * sin_rlat - sin_plat * cos_rlat * cos_rlon
    ) + np.deg2rad(pole_lon)

    return np.rad2deg(lon), np.rad2deg(lat)