# make_ptracer01_ic.py
# Generates PTRACER01 initial condition from real Arctic microplastic sample locations.
# All 68 samples used -- points falling on land will be zeroed by hFacC at runtime.
# Equal weighting, Gaussian spread of SPREAD_DEG around each sample.
# Run from exp2_ptracers/
#
# Output: run/input/ptracers_ic_tr1.bin
#         shape (1, 270, 270), big-endian float64

import numpy as np
import pandas as pd
import netCDF4 as nc

GRID_NC  = "../shared/original_data/ECCO-GRID_06.nc"
CSV      = "../shared/original_data/arctic_mp_samples.csv"
OUT_FILE = "run/input/ptracers_ic_tr1.bin"

# Gaussian spread radius in degrees -- ~100 km at 70N
SPREAD_DEG = 1.0
AMPLITUDE  = 1.0

# Load grid coordinates
with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)  # (270, 270)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)  # (270, 270)

print(f"Grid loaded: shape {XC.shape}")

# Load sample points
df = pd.read_csv(CSV)
print(f"Samples loaded: {len(df)} points")

tracer = np.zeros((270, 270), dtype=np.float64)

for _, row in df.iterrows():
    lat = row["Sample Latitude"]
    lon = row["Sample Longitude"]
    dist2 = (YC - lat)**2 + (XC - lon)**2
    tracer += AMPLITUDE * np.exp(-dist2 / (2.0 * SPREAD_DEG**2))

print(f"Tracer min={tracer.min():.4f}  max={tracer.max():.4f}  nonzero={np.count_nonzero(tracer)}")

# Write big-endian float64, shape (1, 270, 270)
tracer_out = tracer[np.newaxis, :, :].astype('>f8')
tracer_out.tofile(OUT_FILE)

expected = 1 * 270 * 270 * 8
actual   = tracer_out.nbytes
print(f"Wrote {OUT_FILE}  ({actual} bytes, expected {expected})")
