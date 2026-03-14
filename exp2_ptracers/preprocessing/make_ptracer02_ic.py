# make_ptracer02_ic.py
# Generates PTRACER02 initial condition from known Arctic gyre and drift locations.
# Used to verify circulation behavior -- points seeded at known dynamical features.
# Run from exp2_ptracers/
#
# Output: run/input/ptracers_ic_tr2.bin
#         shape (1, 270, 270), big-endian float64

import numpy as np
import netCDF4 as nc

GRID_NC  = "../shared/original_data/ECCO-GRID_06.nc"
OUT_FILE = "run/input/ptracers_ic_tr2.bin"

# Gaussian spread radius in degrees
SPREAD_DEG = 1.0
AMPLITUDE  = 1.0

# Placeholder points -- known Arctic dynamical features
# Format: (lat_deg, lon_deg, label)
SOURCES = [
    (75.0, -140.0, "Beaufort Gyre center"),
    (85.0,   90.0, "Transpolar Drift entry"),
    (79.0,    5.0, "Fram Strait"),
]

# Load grid coordinates
with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)  # (270, 270)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)  # (270, 270)

print(f"Grid loaded: shape {XC.shape}")

tracer = np.zeros((270, 270), dtype=np.float64)

for lat, lon, label in SOURCES:
    dist2 = (YC - lat)**2 + (XC - lon)**2
    tracer += AMPLITUDE * np.exp(-dist2 / (2.0 * SPREAD_DEG**2))
    print(f"  Added: {label} at ({lat}, {lon})")

print(f"Tracer min={tracer.min():.4f}  max={tracer.max():.4f}  nonzero={np.count_nonzero(tracer)}")

# Write big-endian float64, shape (1, 270, 270)
tracer_out = tracer[np.newaxis, :, :].astype('>f8')
tracer_out.tofile(OUT_FILE)

expected = 1 * 270 * 270 * 8
actual   = tracer_out.nbytes
print(f"Wrote {OUT_FILE}  ({actual} bytes, expected {expected})")
