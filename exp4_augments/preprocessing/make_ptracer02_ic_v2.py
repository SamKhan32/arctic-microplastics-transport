# make_ptracer02_ic.py
# Generates PTRACER02 initial condition -- single point at Beaufort Gyre center.
# Snapped to nearest wet cell, value 1.0.
# Run from exp4_wind_fix/
#
# Output: run/input/ptracers_ic_tr2.bin
#         shape (1, 270, 270), big-endian float64

import numpy as np
import netCDF4 as nc

GRID_NC  = "../shared/original_data/ECCO-GRID_06.nc"
BATHY    = "../shared/inputs/bathy_arctic_open.bin"
OUT_FILE = "run/input/ptracers_ic_tr2.bin"

# Beaufort Gyre center
GYRE_LAT = 75.0
GYRE_LON = -140.0

# Load grid coordinates
with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)

print(f"Grid loaded: shape {XC.shape}")

# Load bathy to identify wet cells
bathy = np.fromfile(BATHY, dtype='>f8').reshape(270, 270)
wet   = bathy < 0.0
print(f"Wet cells: {wet.sum()} of {wet.size}")

tracer = np.zeros((270, 270), dtype=np.float64)

dist2 = np.full((270, 270), np.inf)
dist2[wet] = (YC[wet] - GYRE_LAT)**2 + (XC[wet] - GYRE_LON)**2
j, i = np.unravel_index(dist2.argmin(), dist2.shape)
tracer[j, i] = 1.0

print(f"Placed gyre tracer at grid ({i}, {j}), lon={XC[j,i]:.2f}, lat={YC[j,i]:.2f}")
print(f"Tracer min={tracer.min():.4f}  max={tracer.max():.4f}  nonzero={np.count_nonzero(tracer)}")

# Write big-endian float64, shape (1, 270, 270)
tracer_out = tracer[np.newaxis, :, :].astype('>f8')
tracer_out.tofile(OUT_FILE)

expected = 1 * 270 * 270 * 8
actual   = tracer_out.nbytes
print(f"Wrote {OUT_FILE}  ({actual} bytes, expected {expected})")
