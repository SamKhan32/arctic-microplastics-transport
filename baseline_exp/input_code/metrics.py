"""
metrics.py - Generate MITgcm metric binary files from curvilinear grid coordinates.

Computes and writes big-endian float32 binaries for:
  DXG, DYG   - grid-face lengths (u-face, v-face) [m]
  DXC, DYC   - cell-center spacings [m]
  RAC        - cell face areas [m^2]
  angleCosC  - cos(angle) of grid angle at cell centers
  angleSinC  - sin(angle) of grid angle at cell centers

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

  DXG[j,i]   = zonal length of southern face of cell (i,j)  [v-point i to i+1]
  DYG[j,i]   = meridional length of western face of cell (i,j) [u-point j to j+1]
  DXC[j,i]   = zonal distance between T(i-1,j) and T(i,j)
  DYC[j,i]   = meridional distance between T(i,j-1) and T(i,j)
  RAC[j,i]   = area of tracer cell (i,j)

Usage:
    python metrics.py

Requires:
    config.py in the same directory (provides make_grid())
    numpy, scipy
"""

import numpy as np
import os
import sys

# ---------------------------------------------------------------------------
# Import grid from config
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
from config import make_grid, NX, NY

# Earth radius [m]
R_EARTH = 6371000.0

OUT_DIR = "."  # write binaries here; change if needed


# ---------------------------------------------------------------------------
# Haversine great-circle distance (vectorized)
# ---------------------------------------------------------------------------
def haversine(lon1, lat1, lon2, lat2):
    """
    Returns great-circle distance [m] between arrays of (lon,lat) points [deg].
    All inputs broadcast against each other.
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# ---------------------------------------------------------------------------
# Grid-angle computation
# ---------------------------------------------------------------------------
def compute_angle(lon, lat):
    """
    Compute the angle (in radians) of the local i-direction relative to true east
    at each tracer cell center. Uses the forward difference along i.

    angleCosC = cos(angle), angleSinC = sin(angle)

    MITgcm convention: angle is the bearing of the +i direction measured
    clockwise from north, or more practically the angle of the local
    x-axis east-of-north. We compute the azimuth of the vector from T(i,j)
    to T(i+1,j) and subtract 90 degrees to get the angle relative to true east.

    Returns cosA, sinA arrays of shape (NY, NX).
    """
    # Pad longitude/latitude at the right boundary by wrapping (periodic) or
    # extrapolating (closed). For a closed basin we extrapolate.
    lon_r = np.concatenate([lon, 2 * lon[:, -1:] - lon[:, -2:-1]], axis=1)  # (NY, NX+1)
    lat_r = np.concatenate([lat, 2 * lat[:, -1:] - lat[:, -2:-1]], axis=1)

    lon1 = np.radians(lon_r[:, :NX])
    lat1 = np.radians(lat_r[:, :NX])
    lon2 = np.radians(lon_r[:, 1:NX+1])
    lat2 = np.radians(lat_r[:, 1:NX+1])

    # Azimuth of the vector from point 1 to point 2
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azimuth = np.arctan2(x, y)  # bearing clockwise from north

    # Angle east of true east = azimuth - pi/2
    angle = azimuth - np.pi / 2

    return np.cos(angle).astype(np.float64), np.sin(angle).astype(np.float64)


# ---------------------------------------------------------------------------
# Write big-endian float32 binary
# ---------------------------------------------------------------------------
def write_bin(filename, data):
    """Write 2-D array as big-endian float32 binary (MITgcm default)."""
    path = os.path.join(OUT_DIR, filename)
    data.astype(">f4").tofile(path)
    print(f"  Wrote {path}  shape={data.shape}  "
          f"min={data.min():.4g}  max={data.max():.4g}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading grid from config.py ...")
    LON, LAT = make_grid()          # (NY, NX) geographic coordinates [deg]
    assert LON.shape == (NY, NX), f"Expected ({NY},{NX}), got {LON.shape}"
    print(f"  Grid shape: NY={NY}, NX={NX}")

    # ------------------------------------------------------------------
    # Build extended coordinate arrays for differencing.
    # We need values at i+/-0.5 and j+/-0.5 positions, obtained by averaging
    # neighbouring tracer points (corner/node values) and by ghost cells
    # at boundaries (extrapolation for closed basin).
    #
    # Naming: _ext_i has shape (NY, NX+1) - one extra column on the right
    #         _ext_j has shape (NY+1, NX) - one extra row on the top
    # ------------------------------------------------------------------

    # --- i-extended (add ghost column at right) ---
    lon_ei = np.concatenate([LON,
                              2 * LON[:, -1:] - LON[:, -2:-1]], axis=1)  # (NY, NX+1)
    lat_ei = np.concatenate([LAT,
                              2 * LAT[:, -1:] - LAT[:, -2:-1]], axis=1)

    # --- j-extended (add ghost row at top) ---
    lon_ej = np.concatenate([LON,
                              2 * LON[-1:, :] - LON[-2:-1, :]], axis=0)  # (NY+1, NX)
    lat_ej = np.concatenate([LAT,
                              2 * LAT[-1:, :] - LAT[-2:-1, :]], axis=0)

    # ------------------------------------------------------------------
    # DXC[j,i] - zonal distance between T(i-1,j) and T(i,j)
    #   = haversine(lon[j,i-1],lat[j,i-1], lon[j,i],lat[j,i])
    # Left-shifted: prepend ghost column on the left.
    # ------------------------------------------------------------------
    print("Computing DXC ...")
    lon_li = np.concatenate([2 * LON[:, :1] - LON[:, 1:2], LON], axis=1)  # (NY, NX+1)
    lat_li = np.concatenate([2 * LAT[:, :1] - LAT[:, 1:2], LAT], axis=1)
    DXC = haversine(lon_li[:, :NX], lat_li[:, :NX],
                    lon_li[:, 1:],  lat_li[:, 1:])

    # ------------------------------------------------------------------
    # DYC[j,i] - meridional distance between T(i,j-1) and T(i,j)
    # Bottom-shifted: prepend ghost row on the bottom.
    # ------------------------------------------------------------------
    print("Computing DYC ...")
    lon_bj = np.concatenate([2 * LON[:1, :] - LON[1:2, :], LON], axis=0)  # (NY+1, NX)
    lat_bj = np.concatenate([2 * LAT[:1, :] - LAT[1:2, :], LAT], axis=0)
    DYC = haversine(lon_bj[:NY, :], lat_bj[:NY, :],
                    lon_bj[1:,  :], lat_bj[1:,  :])

    # ------------------------------------------------------------------
    # DXG[j,i] - zonal length of the SOUTHERN face of cell (i,j)
    #   = distance between corner points (i-0.5, j-0.5) and (i+0.5, j-0.5)
    # Corner/node lon/lat at (j-0.5, i-0.5) is the average of the four
    # surrounding T points. For the southern face we need nodes along row j-0.5.
    #
    # Node at (j-0.5, i+0.5) ~ average of T(i,j-1), T(i+1,j-1), T(i,j), T(i+1,j)
    # build the node grid first
    # ------------------------------------------------------------------
    print("Computing node (corner) coordinates ...")
    # Extend LON/LAT by one in each direction for node averaging
    # Shape after extension: (NY+1, NX+1)
    lon_x = np.zeros((NY + 1, NX + 1))
    lat_x = np.zeros((NY + 1, NX + 1))

    # Interior nodes: average of 4 surrounding T points
    lon_x[1:NY, 1:NX] = 0.25 * (LON[:NY-1, :NX-1] + LON[:NY-1, 1:NX] +
                                   LON[1:NY,  :NX-1] + LON[1:NY,  1:NX])
    lat_x[1:NY, 1:NX] = 0.25 * (LAT[:NY-1, :NX-1] + LAT[:NY-1, 1:NX] +
                                   LAT[1:NY,  :NX-1] + LAT[1:NY,  1:NX])

    # Boundary nodes - extrapolate
    lon_x[0,  1:NX] = 2 * lon_x[1,  1:NX] - lon_x[2,  1:NX]
    lon_x[NY, 1:NX] = 2 * lon_x[NY-1, 1:NX] - lon_x[NY-2, 1:NX]
    lat_x[0,  1:NX] = 2 * lat_x[1,  1:NX] - lat_x[2,  1:NX]
    lat_x[NY, 1:NX] = 2 * lat_x[NY-1, 1:NX] - lat_x[NY-2, 1:NX]

    lon_x[1:NY, 0]  = 2 * lon_x[1:NY, 1]  - lon_x[1:NY, 2]
    lon_x[1:NY, NX] = 2 * lon_x[1:NY, NX-1] - lon_x[1:NY, NX-2]
    lat_x[1:NY, 0]  = 2 * lat_x[1:NY, 1]  - lat_x[1:NY, 2]
    lat_x[1:NY, NX] = 2 * lat_x[1:NY, NX-1] - lat_x[1:NY, NX-2]

    # Corners
    for jr, jc in [(0, 0), (0, NX), (NY, 0), (NY, NX)]:
        jr2 = 1 if jr == 0 else NY - 1
        jc2 = 1 if jc == 0 else NX - 1
        lon_x[jr, jc] = lon_x[jr2, jc2]
        lat_x[jr, jc] = lat_x[jr2, jc2]

    # DXG[j,i] = distance between node(j, i) and node(j, i+1)  [southern face]
    print("Computing DXG ...")
    DXG = haversine(lon_x[:NY, :NX], lat_x[:NY, :NX],
                    lon_x[:NY, 1:],  lat_x[:NY, 1:])

    # DYG[j,i] = distance between node(j, i) and node(j+1, i)  [western face]
    print("Computing DYG ...")
    DYG = haversine(lon_x[:NY, :NX], lat_x[:NY, :NX],
                    lon_x[1:,  :NX], lat_x[1:,  :NX])

    # RAC[j,i] - area of tracer cell (i,j)
    # Computed as the spherical quadrilateral area from the 4 corner nodes.
    # For efficiency we could use the cross-product approximation:
    #   RAC ~ DXG * DYG  (valid for small cells, ~0.1% error at 10 km)
    # But we use the spherical excess formula for accuracy.
    print("Computing RAC ...")
    # Accurate method: spherical excess of the quadrilateral
    # Split into 2 triangles and sum.
    def spherical_triangle_area(lon_a, lat_a, lon_b, lat_b, lon_c, lat_c):
        """Spherical excess area of triangle ABC [m^2]."""
        a = np.radians(lon_a); A = np.radians(lat_a)
        b = np.radians(lon_b); B = np.radians(lat_b)
        c = np.radians(lon_c); C = np.radians(lat_c)

        # Convert to Cartesian unit vectors
        xa = np.array([np.cos(A)*np.cos(a), np.cos(A)*np.sin(a), np.sin(A)])
        xb = np.array([np.cos(B)*np.cos(b), np.cos(B)*np.sin(b), np.sin(B)])
        xc = np.array([np.cos(C)*np.cos(c), np.cos(C)*np.sin(c), np.sin(C)])

        # Side lengths (angular)
        ab = np.arccos(np.clip(np.sum(xa * xb, axis=0), -1, 1))
        bc = np.arccos(np.clip(np.sum(xb * xc, axis=0), -1, 1))
        ca = np.arccos(np.clip(np.sum(xc * xa, axis=0), -1, 1))

        # Spherical excess via L'Huilier's theorem
        s = 0.5 * (ab + bc + ca)
        tan_e4 = np.sqrt(np.clip(
            np.tan(s/2) * np.tan((s-ab)/2) * np.tan((s-bc)/2) * np.tan((s-ca)/2),
            0, None))
        E = 4 * np.arctan(tan_e4)
        return R_EARTH * R_EARTH * E

    # Quadrilateral corners: SW, SE, NE, NW
    SW_lon, SW_lat = lon_x[:NY, :NX], lat_x[:NY, :NX]
    SE_lon, SE_lat = lon_x[:NY, 1:],  lat_x[:NY, 1:]
    NE_lon, NE_lat = lon_x[1:,  1:],  lat_x[1:,  1:]
    NW_lon, NW_lat = lon_x[1:,  :NX], lat_x[1:,  :NX]

    A1 = spherical_triangle_area(SW_lon, SW_lat, SE_lon, SE_lat, NE_lon, NE_lat)
    A2 = spherical_triangle_area(SW_lon, SW_lat, NE_lon, NE_lat, NW_lon, NW_lat)
    RAC = A1 + A2

    # angleCosC, angleSinC
    print("Computing grid angles ...")
    angleCosC, angleSinC = compute_angle(LON, LAT)

    # ------------------------------------------------------------------
    # Write all files
    # ------------------------------------------------------------------
    print("\nWriting binary files ...")
    write_bin("DXG.bin",       DXG.astype(np.float32))
    write_bin("DYG.bin",       DYG.astype(np.float32))
    write_bin("DXC.bin",       DXC.astype(np.float32))
    write_bin("DYC.bin",       DYC.astype(np.float32))
    write_bin("RAC.bin",       RAC.astype(np.float32))
    write_bin("angleCosC.bin", angleCosC.astype(np.float32))
    write_bin("angleSinC.bin", angleSinC.astype(np.float32))

    # ------------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------------
    print("\nSanity checks:")
    dx_km = DXC.mean() / 1e3
    dy_km = DYC.mean() / 1e3
    area_km2 = RAC.mean() / 1e6
    target_area = 10.7 * 10.7
    print(f"  Mean DXC = {dx_km:.2f} km   (target ~10.7 km)")
    print(f"  Mean DYC = {dy_km:.2f} km   (target ~10.7 km)")
    print(f"  Mean RAC = {area_km2:.1f} km^2  (target ~{target_area:.0f} km^2)")
    angle_err = np.max(np.abs(angleCosC * angleCosC + angleSinC * angleSinC - 1))
    print(f"  cos^2 + sin^2 max deviation: {angle_err:.2e}  (expect ~0)")
    print(f"  RAC negative cells: {(RAC <= 0).sum()}  (expect 0)")
    print(f"  DXG min/max: {DXG.min()/1e3:.2f} / {DXG.max()/1e3:.2f} km")
    print(f"  DYG min/max: {DYG.min()/1e3:.2f} / {DYG.max()/1e3:.2f} km")
    print("\nDone.")


if __name__ == "__main__":
    main()