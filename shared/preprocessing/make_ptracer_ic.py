import numpy as np
import pandas as pd
import netCDF4 as nc
import argparse
import os
import sys

# Usage:
#   tr1 from CSV:   python make_ptracer_ic.py --out run/input/ptracers_ic_tr1.bin --mode csv --nr 5
#   tr2 point src:  python make_ptracer_ic.py --out run/input/ptracers_ic_tr2.bin --mode point --lat 75.0 --lon -140.0 --nr 5

GRID_NC = "../shared/original_data/ECCO-GRID_06.nc"
CSV     = "../shared/original_data/arctic_mp_samples.csv"
BATHY   = "../shared/inputs/bathy_arctic_open.bin"

NY, NX = 270, 270

parser = argparse.ArgumentParser()
parser.add_argument('--out',  required=True, help='output .bin filename')
parser.add_argument('--mode', required=True, choices=['csv', 'point'], help='csv = snap samples, point = single point source')
parser.add_argument('--nr',   type=int, default=1, help='number of vertical levels (default 1)')
parser.add_argument('--lat',  type=float, help='latitude for point mode')
parser.add_argument('--lon',  type=float, help='longitude for point mode')
parser.add_argument('--max-dist', type=float, default=2.0, help='max snap distance in degrees (csv mode, default 2.0)')
args = parser.parse_args()

if os.path.exists(args.out):
    print('ERROR: %s already exists. Use a different --out name.' % args.out)
    sys.exit(1)

if args.mode == 'point' and (args.lat is None or args.lon is None):
    print('ERROR: --lat and --lon required for point mode')
    sys.exit(1)

# Load grid
with nc.Dataset(GRID_NC) as ds:
    XC = np.array(ds.variables['XC'][:], dtype=np.float64)
    YC = np.array(ds.variables['YC'][:], dtype=np.float64)

print('Grid loaded: shape %s' % str(XC.shape))

# Load bathy
bathy = np.fromfile(BATHY, dtype='>f8').reshape(NY, NX)
wet = bathy < 0.0
print('Wet cells: %d of %d' % (wet.sum(), wet.size))

layer = np.zeros((NY, NX), dtype=np.float64)

if args.mode == 'csv':
    df = pd.read_csv(CSV)
    print('Samples loaded: %d points' % len(df))
    placed = 0
    skipped = 0
    for _, row in df.iterrows():
        lat = row['Sample Latitude']
        lon = row['Sample Longitude']
        dist2 = np.full((NY, NX), np.inf)
        dist2[wet] = (YC[wet] - lat)**2 + (XC[wet] - lon)**2
        j, i = np.unravel_index(dist2.argmin(), dist2.shape)
        if dist2[j, i] == np.inf:
            print('  SKIP: no wet cell for (%.2f, %.2f)' % (lat, lon))
            skipped += 1
            continue
        if np.sqrt(dist2[j, i]) > args.max_dist:
            print('  SKIP: (%.2f, %.2f) snaps to (%.2f, %.2f), too far' % (lat, lon, XC[j,i], YC[j,i]))
            skipped += 1
            continue
        layer[j, i] = 1.0
        placed += 1
    print('Placed %d samples, skipped %d' % (placed, skipped))

elif args.mode == 'point':
    dist2 = np.full((NY, NX), np.inf)
    dist2[wet] = (YC[wet] - args.lat)**2 + (XC[wet] - args.lon)**2
    j, i = np.unravel_index(dist2.argmin(), dist2.shape)
    if dist2[j, i] == np.inf:
        print('ERROR: no wet cell found near (%.2f, %.2f)' % (args.lat, args.lon))
        sys.exit(1)
    layer[j, i] = 1.0
    print('Point placed at grid (j=%d, i=%d) -> (%.2f, %.2f)' % (j, i, YC[j,i], XC[j,i]))

# Build 3D array: tracer in level 0 only, zeros below
out = np.zeros((args.nr, NY, NX), dtype='>f8')
out[0, :, :] = layer

print('Tracer nonzero cells: %d  max=%.4f' % (np.count_nonzero(layer), layer.max()))

out.tofile(args.out)
expected = args.nr * NY * NX * 8
print('Wrote %s  (%d bytes, expected %d)' % (args.out, out.nbytes, expected))
