"""
windStress.py - Download ERA5 monthly mean surface stress, compute a 12-month
climatology, interpolate onto the curvilinear grid, rotate into grid-relative
coordinates, and write MITgcm-ready big-endian float32 binaries.

Pipeline:
    1. Download metss/mntss from CDS (1991-2020, all months) -> original_data/
    2. Average over years to get climatological monthly means (12 x NY x NX)
    3. Interpolate from ERA5 0.25-deg grid onto curvilinear LON/LAT
    4. Rotate (tau_east, tau_north) into grid-relative (tau_x, tau_y)
       using angleCosC/angleSinC from metrics.py output
    5. Write oceTauX_mm.bin, oceTauY_mm.bin for mm = 01 .. 12

Output files (in MITGCM_INPUT_DIR):
    oceTauX_01.bin ... oceTauX_12.bin   [N/m^2, grid-relative x]
    oceTauY_01.bin ... oceTauY_12.bin   [N/m^2, grid-relative y]

Usage:
    python windStress.py           # download + process
    python windStress.py --no-download  # skip download, reuse existing NetCDF
"""

import argparse
import os
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import make_grid, NX, NY

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT_DIR        = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ORIGINAL_DIR    = os.path.join(ROOT_DIR, "original_data")
MITGCM_DIR      = os.path.join(ROOT_DIR, "mitgcm_input")
ERA5_FILE       = os.path.join(ORIGINAL_DIR, "era5_stress_monthly_1991_2020.nc")

os.makedirs(ORIGINAL_DIR, exist_ok=True)
os.makedirs(MITGCM_DIR,   exist_ok=True)

# ---------------------------------------------------------------------------
# CDS download
# ---------------------------------------------------------------------------
def download_era5():
    print("Connecting to CDS ...")
    import cdsapi
    c = cdsapi.Client()

    years  = [str(y) for y in range(1991, 2021)]
    months = [f"{m:02d}" for m in range(1, 13)]

    print(f"Requesting ERA5 surface stress 1991-2020 ({len(years)*len(months)} months) ...")
    print("This may queue for a few minutes on the CDS side.")

    c.retrieve(
        "reanalysis-era5-single-levels-monthly-means",
        {
            "product_type": "monthly_averaged_reanalysis",
            "variable": [
                "mean_eastward_turbulent_surface_stress",
                "mean_northward_turbulent_surface_stress",
            ],
            "year":   years,
            "month":  months,
            "time":   "00:00",
            "format": "netcdf",
        },
        ERA5_FILE,
    )
    print(f"Downloaded -> {ERA5_FILE}")


# ---------------------------------------------------------------------------
# Load and build climatology
# ---------------------------------------------------------------------------
def load_climatology():
    """
    Read ERA5 NetCDF, average over years, return:
        era_lon  : (nlon,)       0..360 or -180..180
        era_lat  : (nlat,)       descending (90..-90) as ERA5 delivers
        tau_e_clim : (12, nlat, nlon)  eastward stress [N/m^2]
        tau_n_clim : (12, nlat, nlon)  northward stress [N/m^2]
    """
    try:
        import netCDF4 as nc
    except ImportError:
        raise ImportError("netCDF4 not installed. Run: pip install netCDF4")

    print(f"Loading {ERA5_FILE} ...")
    ds = nc.Dataset(ERA5_FILE)

    # ERA5 variable names in the monthly-means product
    # metss = mean eastward turbulent surface stress
    # mntss = mean northward turbulent surface stress
    tau_e_raw = ds.variables["avg_iews"][:]   # (time, lat, lon)
    tau_n_raw = ds.variables["avg_inss"][:]
    era_lat   = ds.variables["latitude"][:]
    era_lon   = ds.variables["longitude"][:]
    ds.close()

    # time axis is (year*12 + month) ordered, 360 entries for 30 years
    ntimes = tau_e_raw.shape[0]
    assert ntimes == 360, f"Expected 360 time steps, got {ntimes}"

    # reshape to (nyear, 12, nlat, nlon) then mean over years
    tau_e = tau_e_raw.reshape(30, 12, len(era_lat), len(era_lon))
    tau_n = tau_n_raw.reshape(30, 12, len(era_lat), len(era_lon))

    tau_e_clim = tau_e.mean(axis=0)   # (12, nlat, nlon)
    tau_n_clim = tau_n.mean(axis=0)
    print(f"  ERA5 raw tau_e abs max: {np.abs(tau_e_clim).max():.4f}")

    print(f"  ERA5 grid: {len(era_lat)} lat x {len(era_lon)} lon")
    print(f"  Climatology range tau_e: {tau_e_clim.min():.4f} to {tau_e_clim.max():.4f} N/m^2")
    print(f"  Climatology range tau_n: {tau_n_clim.min():.4f} to {tau_n_clim.max():.4f} N/m^2")

    return np.array(era_lon), np.array(era_lat), tau_e_clim, tau_n_clim


# ---------------------------------------------------------------------------
# Interpolation onto curvilinear grid
# ---------------------------------------------------------------------------
def interpolate_to_grid(era_lon, era_lat, tau_e_clim, tau_n_clim, LON, LAT):
    """
    Bilinear interpolation from ERA5 regular grid onto curvilinear (LON, LAT).

    ERA5 lon is 0..359.75. We need to handle the wrap at 360/0 for points
    in the curvilinear grid that may be expressed as negative longitudes.

    Returns tau_e, tau_n each of shape (12, NY, NX).
    """
    print("Interpolating onto curvilinear grid ...")

    # Normalise target longitudes to 0..360 to match ERA5
    LON_360 = LON % 360.0

    # ERA5 lat is descending (90 -> -90), RegularGridInterpolator needs ascending
    # so we flip and note the index order
    lat_asc = era_lat[::-1]   # ascending
    nlat = len(lat_asc)
    nlon = len(era_lon)

    # Extend ERA5 longitude by one point on each side to handle wrap-around
    # during interpolation near the 0/360 boundary
    lon_ext = np.concatenate([era_lon[-1:] - 360.0, era_lon, era_lon[:1] + 360.0])

    tau_e_out = np.zeros((12, NY, NX), dtype=np.float32)
    tau_n_out = np.zeros((12, NY, NX), dtype=np.float32)

    for month in range(12):
        # Flip lat axis to ascending
        field_e = tau_e_clim[month, ::-1, :]   # (nlat, nlon)
        field_n = tau_n_clim[month, ::-1, :]

        # Extend longitude dimension for wrap
        field_e_ext = np.concatenate([field_e[:, -1:], field_e, field_e[:, :1]], axis=1)
        field_n_ext = np.concatenate([field_n[:, -1:], field_n, field_n[:, :1]], axis=1)

        interp_e = RegularGridInterpolator(
            (lat_asc, lon_ext), field_e_ext,
            method="linear", bounds_error=False, fill_value=None
        )
        interp_n = RegularGridInterpolator(
            (lat_asc, lon_ext), field_n_ext,
            method="linear", bounds_error=False, fill_value=None
        )

        pts = np.column_stack([LAT.ravel(), LON_360.ravel()])
        tau_e_out[month] = interp_e(pts).reshape(NY, NX)
        tau_n_out[month] = interp_n(pts).reshape(NY, NX)

        print(f"  Month {month+1:02d} done", end="\r")

    print()
    return tau_e_out, tau_n_out


# ---------------------------------------------------------------------------
# Rotate into grid-relative coordinates
# ---------------------------------------------------------------------------
def rotate_to_grid(tau_e, tau_n, cos_angle, sin_angle):
    """
    Rotate geographic (east, north) stress components into grid-relative (x, y).

    The rotation uses angleCosC and angleSinC from metrics.py:
        tau_x =  tau_e * cos_angle + tau_n * sin_angle
        tau_y = -tau_e * sin_angle + tau_n * cos_angle

    This is the standard rotation for MITgcm curvilinear grids.
    tau_x drives the u-momentum equation, tau_y drives v-momentum.

    Args:
        tau_e, tau_n  : (12, NY, NX) geographic stress [N/m^2]
        cos_angle     : (NY, NX) angleCosC from metrics output
        sin_angle     : (NY, NX) angleSinC from metrics output

    Returns:
        tau_x, tau_y  : (12, NY, NX) grid-relative stress [N/m^2]
    """
    print("Rotating stress into grid-relative coordinates ...")
    tau_x =  tau_e * cos_angle[np.newaxis, :, :] + tau_n * sin_angle[np.newaxis, :, :]
    tau_y = -tau_e * sin_angle[np.newaxis, :, :] + tau_n * cos_angle[np.newaxis, :, :]
    return tau_x, tau_y


# ---------------------------------------------------------------------------
# Load angle files from metrics output
# ---------------------------------------------------------------------------
def load_angles():
    cos_path = os.path.join(MITGCM_DIR, "angleCosC.bin")
    sin_path = os.path.join(MITGCM_DIR, "angleSinC.bin")
    if not os.path.exists(cos_path) or not os.path.exists(sin_path):
        raise FileNotFoundError(
            "angleCosC.bin / angleSinC.bin not found in mitgcm_input/. "
            "Run metrics.py first."
        )
    cos_angle = np.fromfile(cos_path, dtype=">f4").reshape(NY, NX)
    sin_angle = np.fromfile(sin_path, dtype=">f4").reshape(NY, NX)
    return cos_angle, sin_angle


# ---------------------------------------------------------------------------
# Write binary files
# ---------------------------------------------------------------------------
def write_bin(filename, data):
    path = os.path.join(MITGCM_DIR, filename)
    data.astype(">f4").tofile(path)
    print(f"  Wrote {path}  min={data.min():.4f}  max={data.max():.4f} N/m^2")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-download", action="store_true",
        help="Skip CDS download and use existing NetCDF file"
    )
    args = parser.parse_args()

    # Step 1: download
    if args.no_download:
        if not os.path.exists(ERA5_FILE):
            raise FileNotFoundError(f"--no-download set but file not found: {ERA5_FILE}")
        print(f"Skipping download, using {ERA5_FILE}")
    else:
        download_era5()

    # Step 2: climatology
    era_lon, era_lat, tau_e_clim, tau_n_clim = load_climatology()

    # Step 3: load grid
    print("Loading curvilinear grid ...")
    LON, LAT = make_grid()

    # Step 4: interpolate
    tau_e, tau_n = interpolate_to_grid(era_lon, era_lat, tau_e_clim, tau_n_clim, LON, LAT)

    # Step 5: load angles and rotate
    cos_angle, sin_angle = load_angles()
    tau_x, tau_y = rotate_to_grid(tau_e, tau_n, cos_angle, sin_angle)
    # mask wind stress to ocean points only
    bathy_path = os.path.join(MITGCM_DIR, "bathy.bin")
    bathy = np.fromfile(bathy_path, dtype=">f4").reshape(NY, NX)
    ocean_mask = bathy != 0.0   # True where ocean

    for month in range(12):
        tau_x[month] *= ocean_mask
        tau_y[month] *= ocean_mask
    # Step 6: write
    print("\nWriting binary files ...")
    tau_x.astype(">f4").tofile(os.path.join(MITGCM_DIR, "oceTauX.bin"))
    tau_y.astype(">f4").tofile(os.path.join(MITGCM_DIR, "oceTauY.bin"))


    # Sanity checks
    print("\nSanity checks:")
    print(f"  tau_x annual mean: {tau_x.mean():.4f} N/m^2  (expect small, near zero)")
    print(f"  tau_y annual mean: {tau_y.mean():.4f} N/m^2  (expect small, near zero)")
    print(f"  tau_x max magnitude: {np.abs(tau_x).max():.4f} N/m^2  (expect < 1.0)")
    print(f"  tau_y max magnitude: {np.abs(tau_y).max():.4f} N/m^2  (expect < 1.0)")
    print("\nDone.")


if __name__ == "__main__":
    main()
