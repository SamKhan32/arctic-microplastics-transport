"""
metrics.py - Generate MITgcm metric binary files from curvilinear grid coordinates.

Computes all fields required by horizGridFile (18 records) plus writes them
individually for reference. Record order matches ini_curvilinear_grid.F:

  Record  1  : xC        - tracer cell center longitude [deg]
  Record  2  : yC        - tracer cell center latitude [deg]
  Record  3  : dxF       - cell face length in x at cell center [m]
  Record  4  : dyF       - cell face length in y at cell center [m]
  Record  5  : rA        - tracer cell area [m^2]
  Record  6  : xG        - corner (vorticity) longitude [deg]
  Record  7  : yG        - corner (vorticity) latitude [deg]
  Record  8  : dxV       - x face length at vorticity points [m]
  Record  9  : dyU       - y face length at u points [m]
  Record  10 : rAz       - vorticity cell area [m^2]
  Record  11 : dxC       - x distance between tracer cell centers [m]
  Record  12 : dyC       - y distance between tracer cell centers [m]
  Record  13 : rAw       - u-face cell area [m^2]
  Record  14 : rAs       - v-face cell area [m^2]
  Record  15 : dxG       - x length of southern cell face [m]
  Record  16 : dyG       - y length of western cell face [m]
  Record  17 : angleCosC - cos of grid angle at cell centers
  Record  18 : angleSinC - sin of grid angle at cell centers

Grid staggering (MITgcm Arakawa C):

        v(i,j+1)
         |
  u(i,j)--T(i,j)--u(i+1,j)
         |
        v(i,j)

  T (tracer)   : cell centers      - LON[j,i], LAT[j,i]
  u            : west face centers  - between T(i-1,j) and T(i,j)
  v            : south face centers - between T(i,j-1) and T(i,j)
  vorticity    : corner/node points - between 4 surrounding T cells

Usage:
    python metrics.py

Requires:
    config.py in repo root (provides make_grid())
    numpy
"""

from pathlib import Path
import numpy as np
import os
import sys

# ---------------------------------------------------------------------------
# Import grid from config
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import make_grid, NX, NY

R_EARTH = 6371000.0

ROOT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR  = str(ROOT_DIR / "mitgcm_input")


# ---------------------------------------------------------------------------
# Haversine great-circle distance (vectorized)
# ---------------------------------------------------------------------------
def haversine(lon1, lat1, lon2, lat2):
    """Great-circle distance [m] between arrays of (lon, lat) points [deg]."""
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (np.sin(dlat / 2) * np.sin(dlat / 2)
         + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) * np.sin(dlon / 2))
    return 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# ---------------------------------------------------------------------------
# Spherical quadrilateral area via L'Huilier's theorem (two triangles)
# ---------------------------------------------------------------------------
def spherical_triangle_area(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c):
    """Spherical excess area of triangle ABC [m^2]."""
    a = np.radians(lon_a); A = np.radians(lat_a)
    b = np.radians(lon_b); B = np.radians(lat_b)
    c = np.radians(lon_c); C = np.radians(lat_c)

    xa = np.array([np.cos(A)*np.cos(a), np.cos(A)*np.sin(a), np.sin(A)])
    xb = np.array([np.cos(B)*np.cos(b), np.cos(B)*np.sin(b), np.sin(B)])
    xc = np.array([np.cos(C)*np.cos(c), np.cos(C)*np.sin(c), np.sin(C)])

    ab = np.arccos(np.clip(np.sum(xa * xb, axis=0), -1, 1))
    bc = np.arccos(np.clip(np.sum(xb * xc, axis=0), -1, 1))
    ca = np.arccos(np.clip(np.sum(xc * xa, axis=0), -1, 1))

    s = 0.5 * (ab + bc + ca)
    tan_e4 = np.sqrt(np.clip(
        np.tan(s/2) * np.tan((s-ab)/2) * np.tan((s-bc)/2) * np.tan((s-ca)/2),
        0, None))
    E = 4 * np.arctan(tan_e4)
    return R_EARTH * R_EARTH * E


def quad_area(lon_sw, lat_sw, lon_se, lat_se, lon_ne, lat_ne, lon_nw, lat_nw):
    """Spherical quadrilateral area [m^2] split into two triangles."""
    A1 = spherical_triangle_area(lon_sw, lat_sw, lon_se, lat_se, lon_ne, lat_ne)
    A2 = spherical_triangle_area(lon_sw, lat_sw, lon_ne, lat_ne, lon_nw, lat_nw)
    return A1 + A2


# ---------------------------------------------------------------------------
# Grid-angle computation
# ---------------------------------------------------------------------------
def compute_angle(lon, lat):
    """
    cos and sin of the angle between the local i-direction and true east,
    at each tracer cell center. Shape (NY, NX).
    """
    lon_r = np.concatenate([lon, 2 * lon[:, -1:] - lon[:, -2:-1]], axis=1)
    lat_r = np.concatenate([lat, 2 * lat[:, -1:] - lat[:, -2:-1]], axis=1)

    lon1 = np.radians(lon_r[:, :NX])
    lat1 = np.radians(lat_r[:, :NX])
    lon2 = np.radians(lon_r[:, 1:NX+1])
    lat2 = np.radians(lat_r[:, 1:NX+1])

    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azimuth = np.arctan2(x, y)

    angle = azimuth - np.pi / 2
    return np.cos(angle), np.sin(angle)


# ---------------------------------------------------------------------------
# Write big-endian float64 binary (MITgcm precFloat64)
# ---------------------------------------------------------------------------
def write_bin(filename, data):
    path = os.path.join(OUT_DIR, filename)
    data.astype(">f8").tofile(path)
    print(f"  Wrote {filename:20s}  min={data.min():.4g}  max={data.max():.4g}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading grid from config.py ...")
    LON, LAT = make_grid()
    assert LON.shape == (NY, NX), f"Expected ({NY},{NX}), got {LON.shape}"
    print(f"  Grid shape: NY={NY}, NX={NX}")

    # ------------------------------------------------------------------
    # Corner (vorticity/node) coordinates - shape (NY+1, NX+1)
    # Interior nodes: average of 4 surrounding T points.
    # Boundary nodes: linear extrapolation (closed basin).
    # ------------------------------------------------------------------
    print("Computing corner coordinates ...")
    lon_x = np.zeros((NY + 1, NX + 1))
    lat_x = np.zeros((NY + 1, NX + 1))

    lon_x[1:NY, 1:NX] = 0.25 * (LON[:NY-1, :NX-1] + LON[:NY-1, 1:NX] +
                                   LON[1:NY,  :NX-1] + LON[1:NY,  1:NX])
    lat_x[1:NY, 1:NX] = 0.25 * (LAT[:NY-1, :NX-1] + LAT[:NY-1, 1:NX] +
                                   LAT[1:NY,  :NX-1] + LAT[1:NY,  1:NX])

    lon_x[0,  1:NX] = 2 * lon_x[1,    1:NX] - lon_x[2,    1:NX]
    lon_x[NY, 1:NX] = 2 * lon_x[NY-1, 1:NX] - lon_x[NY-2, 1:NX]
    lat_x[0,  1:NX] = 2 * lat_x[1,    1:NX] - lat_x[2,    1:NX]
    lat_x[NY, 1:NX] = 2 * lat_x[NY-1, 1:NX] - lat_x[NY-2, 1:NX]

    lon_x[1:NY, 0]  = 2 * lon_x[1:NY, 1]    - lon_x[1:NY, 2]
    lon_x[1:NY, NX] = 2 * lon_x[1:NY, NX-1] - lon_x[1:NY, NX-2]
    lat_x[1:NY, 0]  = 2 * lat_x[1:NY, 1]    - lat_x[1:NY, 2]
    lat_x[1:NY, NX] = 2 * lat_x[1:NY, NX-1] - lat_x[1:NY, NX-2]

    for jr, jc in [(0, 0), (0, NX), (NY, 0), (NY, NX)]:
        jr2 = 1 if jr == 0 else NY - 1
        jc2 = 1 if jc == 0 else NX - 1
        lon_x[jr, jc] = lon_x[jr2, jc2]
        lat_x[jr, jc] = lat_x[jr2, jc2]

    # SW corner of cell (j,i) is node (j,i)
    SW_lon = lon_x[:NY, :NX];  SW_lat = lat_x[:NY, :NX]
    SE_lon = lon_x[:NY, 1:];   SE_lat = lat_x[:NY, 1:]
    NE_lon = lon_x[1:,  1:];   NE_lat = lat_x[1:,  1:]
    NW_lon = lon_x[1:,  :NX];  NW_lat = lat_x[1:,  :NX]

    # ------------------------------------------------------------------
    # xC, yC (records 1, 2)
    # xG, yG (records 6, 7)
    # ------------------------------------------------------------------
    xC = LON
    yC = LAT
    xG = SW_lon
    yG = SW_lat

    # ------------------------------------------------------------------
    # dxG, dyG (records 15, 16)
    # dxG = x length of southern face = SW to SE
    # dyG = y length of western face  = SW to NW
    # ------------------------------------------------------------------
    print("Computing dxG, dyG ...")
    dxG = haversine(SW_lon, SW_lat, SE_lon, SE_lat)
    dyG = haversine(SW_lon, SW_lat, NW_lon, NW_lat)

    # ------------------------------------------------------------------
    # dxC, dyC (records 11, 12)
    # dxC = x distance between T(j,i-1) and T(j,i)
    # dyC = y distance between T(j-1,i) and T(j,i)
    # ------------------------------------------------------------------
    print("Computing dxC, dyC ...")
    lon_li = np.concatenate([2 * LON[:, :1] - LON[:, 1:2], LON], axis=1)
    lat_li = np.concatenate([2 * LAT[:, :1] - LAT[:, 1:2], LAT], axis=1)
    dxC = haversine(lon_li[:, :NX], lat_li[:, :NX], lon_li[:, 1:], lat_li[:, 1:]) 

    lon_bj = np.concatenate([2 * LON[:1, :] - LON[1:2, :], LON], axis=0)
    lat_bj = np.concatenate([2 * LAT[:1, :] - LAT[1:2, :], LAT], axis=0)
    dyC = haversine(lon_bj[:NY, :], lat_bj[:NY, :], lon_bj[1:, :], lat_bj[1:, :])
    # After computing dxC, fix northern boundary  
    dxC[-1, :] = dxC[-2, :]
    dxC[:, 0] = dxC[:, 1]
    # ------------------------------------------------------------------
    # dxF, dyF (records 3, 4)
    # dxF = average of northern and southern faces
    # dyF = average of eastern and western faces
    # ------------------------------------------------------------------
    print("Computing dxF, dyF ...")
    dxG_ext = np.concatenate([dxG, 2*dxG[-1:, :] - dxG[-2:-1, :]], axis=0)
    dxF = 0.5 * (dxG_ext[:NY, :] + dxG_ext[1:, :])

    dyG_ext = np.concatenate([dyG, 2*dyG[:, -1:] - dyG[:, -2:-1]], axis=1)
    dyF = 0.5 * (dyG_ext[:, :NX] + dyG_ext[:, 1:])
    # After computing dxF, fix northern boundary
    dxF[-1, :] = dxF[-2, :]

    # ------------------------------------------------------------------
    # dxV, dyU (records 8, 9)
    # dxV = x face length at vorticity points
    # dyU = y face length at u points
    # ------------------------------------------------------------------
    print("Computing dxV, dyU ...")
    dxG_extb = np.concatenate([2*dxG[:1, :] - dxG[1:2, :], dxG], axis=0)
    dxV = 0.5 * (dxG_extb[:NY, :] + dxG_extb[1:, :])

    dyG_extb = np.concatenate([2*dyG[:, :1] - dyG[:, 1:2], dyG], axis=1)
    dyU = 0.5 * (dyG_extb[:, :NX] + dyG_extb[:, 1:])
    # After computing dyU, fix any negatives
    dyU[dyU < 0] = np.abs(dyU[dyU < 0])

    # ------------------------------------------------------------------
    # Cell areas (records 5, 10, 13, 14)
    # rA  = tracer cell area
    # rAz = vorticity cell area
    # rAw = u-face area (approximated as rA at 10km resolution)
    # rAs = v-face area (approximated as rA at 10km resolution)
    # ------------------------------------------------------------------
    print("Computing cell areas ...")

    rA = quad_area(SW_lon, SW_lat, SE_lon, SE_lat, NE_lon, NE_lat, NW_lon, NW_lat)

    # rAz: vorticity cell at corner(j,i), bounded by T(j-1,i-1), T(j-1,i), T(j,i), T(j,i-1)
    lon_ext = np.zeros((NY+2, NX+2))
    lat_ext = np.zeros((NY+2, NX+2))
    lon_ext[1:NY+1, 1:NX+1] = LON
    lat_ext[1:NY+1, 1:NX+1] = LAT
    lon_ext[0,  1:NX+1] = 2*LON[0,  :] - LON[1,  :]
    lon_ext[-1, 1:NX+1] = 2*LON[-1, :] - LON[-2, :]
    lat_ext[0,  1:NX+1] = 2*LAT[0,  :] - LAT[1,  :]
    lat_ext[-1, 1:NX+1] = 2*LAT[-1, :] - LAT[-2, :]
    lon_ext[1:NY+1, 0]  = 2*LON[:, 0]  - LON[:, 1]
    lon_ext[1:NY+1, -1] = 2*LON[:, -1] - LON[:, -2]
    lat_ext[1:NY+1, 0]  = 2*LAT[:, 0]  - LAT[:, 1]
    lat_ext[1:NY+1, -1] = 2*LAT[:, -1] - LAT[:, -2]
    lon_ext[0,  0]  = lon_ext[1,  1];   lat_ext[0,  0]  = lat_ext[1,  1]
    lon_ext[0,  -1] = lon_ext[1,  -2];  lat_ext[0,  -1] = lat_ext[1,  -2]
    lon_ext[-1, 0]  = lon_ext[-2, 1];   lat_ext[-1, 0]  = lat_ext[-2, 1]
    lon_ext[-1, -1] = lon_ext[-2, -2];  lat_ext[-1, -1] = lat_ext[-2, -2]

    rAz = quad_area(
        lon_ext[:NY,    :NX],    lat_ext[:NY,    :NX],
        lon_ext[:NY,    1:NX+1], lat_ext[:NY,    1:NX+1],
        lon_ext[1:NY+1, 1:NX+1], lat_ext[1:NY+1, 1:NX+1],
        lon_ext[1:NY+1, :NX],   lat_ext[1:NY+1, :NX],
    )

    rAw = rA.copy()
    rAs = rA.copy()

    # ------------------------------------------------------------------
    # angleCosC, angleSinC (records 17, 18)
    # ------------------------------------------------------------------
    print("Computing grid angles ...")
    angleCosC, angleSinC = compute_angle(LON, LAT)
    print("dxC row 260 (middle), cols 0-5:", dxC[260, :5]/1000)
    print("dxC row 260, cols 764-767:", dxC[260, 764:768]/1000)
    print("LON row 260, cols 0-5:", LON[260, :5])
    print("lon_li row 260, cols 0-5:", lon_li[260, :6])
    print("haversine check:", haversine(lon_li[260:261, :5], lat_li[260:261, :5], 
                                        lon_li[260:261, 1:6], lat_li[260:261, 1:6])/1000)
    # ------------------------------------------------------------------
    # Write individual files (for windStress.py and debugging)
    # ------------------------------------------------------------------
    print("\nWriting individual binary files ...")
    write_bin("DXG.bin",       dxG)
    write_bin("DYG.bin",       dyG)
    write_bin("DXC.bin",       dxC)
    write_bin("DYC.bin",       dyC)
    write_bin("RAC.bin",       rA)
    write_bin("angleCosC.bin", angleCosC)
    write_bin("angleSinC.bin", angleSinC)
    MAX_DX = 50000.0

    for arr in [dxF, dxV, dxG]:
        bad = arr > MAX_DX
        print(f"Bad cells: {bad.sum()}")
        rows, cols = np.where(bad)
        for r, c in zip(rows, cols):
            # average from row above and below, same column
            neighbors = []
            if r > 0 and arr[r-1, c] <= MAX_DX:
                neighbors.append(arr[r-1, c])
            if r < NY-1 and arr[r+1, c] <= MAX_DX:
                neighbors.append(arr[r+1, c])
            if neighbors:
                arr[r, c] = np.mean(neighbors)
            else:
                arr[r, c] = 10700.0  # fallback to ~mean grid spacing
    # ------------------------------------------------------------------
    # Pack tile001.mitgrid in the exact order MITgcm expects (precFloat64)
    # ------------------------------------------------------------------
    print("\nPacking tile001.mitgrid ...")
    horiz_path = os.path.join(OUT_DIR, "tile001.mitgrid")
    records = [
        ( 1, "xC",        xC),
        ( 2, "yC",        yC),
        ( 3, "dxF",       dxF),
        ( 4, "dyF",       dyF),
        ( 5, "rA",        rA),
        ( 6, "xG",        xG),
        ( 7, "yG",        yG),
        ( 8, "dxV",       dxV),
        ( 9, "dyU",       dyU),
        (10, "rAz",       rAz),
        (11, "dxC",       dxC),
        (12, "dyC",       dyC),
        (13, "rAw",       rAw),
        (14, "rAs",       rAs),
        (15, "dxG",       dxG),
        (16, "dyG",       dyG),
        (17, "angleCosC", angleCosC),
        (18, "angleSinC", angleSinC),
    ]
    with open(horiz_path, "wb") as f:
        for rec, name, arr in records:
            f.write(arr.astype(">f8").tobytes())
            print(f"  Record {rec:2d}  {name:12s}  min={arr.min():.4g}  max={arr.max():.4g}")
    print(f"  -> {horiz_path}")

    # ------------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------------
    print("\nSanity checks:")
    print(f"  Mean dxC = {dxC.mean()/1e3:.2f} km   (target ~10.7 km)")
    print(f"  Mean dyC = {dyC.mean()/1e3:.2f} km   (target ~10.7 km)")
    print(f"  Mean rA  = {rA.mean()/1e6:.1f} km^2  (target ~{10.7*10.7:.0f} km^2)")
    angle_err = np.max(np.abs(angleCosC * angleCosC + angleSinC * angleSinC - 1))
    print(f"  cos^2 + sin^2 max deviation: {angle_err:.2e}  (expect ~0)")
    print(f"  rA  negative cells: {(rA  <= 0).sum()}  (expect 0)")
    print(f"  rAz negative cells: {(rAz <= 0).sum()}  (expect 0)")
    total_size_mb = 18 * NY * NX * 8 / 1e6
    print(f"  tile001.mitgrid size: {total_size_mb:.1f} MB")
    print("\nDone.")


if __name__ == "__main__":
    main()