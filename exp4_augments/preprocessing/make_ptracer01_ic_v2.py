# make_ptracer01_ic.py
# Generates PTRACER01 initial condition from real Arctic microplastic sample locations.
# Each sample snapped to nearest wet cell, value 1.0 at that cell only.
# Run from exp4_wind_fix/
#
# Output: run/input/ptracers_ic_tr1.bin
#         shape (1, 270, 270), big-endian float64

import numpy as np
import pandas as pd
import netCDF4 as nc

GRID_NC  = "../shared/original_data/ECCO-GRID_06.nc"
CSV      = "../shared/original_data/arctic_mp_samples.csv"
BATHY    = "../shared/inputs/bathy_arctic_open.bin"
OUT_FILE = "run/input/ptracers_ic_tr1.bin"

# Load grid coordinates
with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)  # (270, 270)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)  # (270, 270)

print(f"Grid loaded: shape {XC.shape}")

# Load bathy to identify wet cells (depth < 0 = ocean)
bathy = np.fromfile(BATHY, dtype='>f8').reshape(270, 270)
wet   = bathy < 0.0
print(f"Wet cells: {wet.sum()} of {wet.size}")

# Load sample points
df = pd.read_csv(CSV)
print(f"Samples loaded: {len(df)} points")

tracer = np.zeros((270, 270), dtype=np.float64)
placed = 0
skipped = 0

for _, row in df.iterrows():
    lat = row["Sample Latitude"]
    lon = row["Sample Longitude"]

    # distance to all wet cells only
    dist2 = np.full((270, 270), np.inf)
    dist2[wet] = (YC[wet] - lat)**2 + (XC[wet] - lon)**2

    j, i = np.unravel_index(dist2.argmin(), dist2.shape)

    if dist2[j, i] == np.inf:
        print(f"  SKIP: no wet cell found for ({lat}, {lon})")
        skipped += 1
        continue

    tracer[j, i] = 1.0
    placed += 1

print(f"Placed {placed} samples, skipped {skipped}")
print(f"Tracer min={tracer.min():.4f}  max={tracer.max():.4f}  nonzero={np.count_nonzero(tracer)}")

# Write big-endian float64, shape (1, 270, 270)
tracer_out = tracer[np.newaxis, :, :].astype('>f8')
tracer_out.tofile(OUT_FILE)

expected = 1 * 270 * 270 * 8
actual   = tracer_out.nbytes
print(f"Wrote {OUT_FILE}  ({actual} bytes, expected {expected})")
