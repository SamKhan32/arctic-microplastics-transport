import netCDF4 as nc
import numpy as np
import os

GRID_PATH = 'original_data/ECCO-GRID_06.nc'
OUT_DIR   = 'mitgcm_input'

Nx,  Ny  = 270, 270
sNx, sNy = 27, 30
nPx, nPy = 10, 9
N_FIELDS  = 18

os.makedirs(OUT_DIR, exist_ok=True)

ds = nc.Dataset(GRID_PATH)

def get(varname):
    return np.array(ds.variables[varname][:], dtype=np.float64)

# 18 fields in the order MITgcm expects in a .mitgrid file
fields_full = np.stack([
    get('XC'),   # 1  xC
    get('YC'),   # 2  yC
    get('dxG'),  # 3  dxF
    get('dyG'),  # 4  dyF
    get('rA'),   # 5  rA
    get('XG'),   # 6  xG
    get('YG'),   # 7  yG
    get('dxC'),  # 8  dxV
    get('dyC'),  # 9  dyU
    get('rAz'),  # 10 rAz
    get('dxC'),  # 11 dxC
    get('dyC'),  # 12 dyC
    get('rAw'),  # 13 rAw
    get('rAs'),  # 14 rAs
    get('dxG'),  # 15 dxG
    get('dyG'),  # 16 dyG
    get('CS'),   # 17 angleCosC
    get('SN'),   # 18 angleSinC
])  # (18, 270, 270)

ds.close()

rank = 0
for j in range(nPy):
    for i in range(nPx):
        rank += 1
        x0 = i * sNx
        y0 = j * sNy
        x1e = min(x0 + sNx + 1, Nx)
        y1e = min(y0 + sNy + 1, Ny)

        tile = np.zeros((N_FIELDS, sNy + 1, sNx + 1), dtype='>f8')
        tile[:, :y1e - y0, :x1e - x0] = fields_full[:, y0:y1e, x0:x1e]

        fname = f'{OUT_DIR}/arctic_cap_v2.mitgrid.face{rank:03d}.bin'
        tile.astype('>f8').tofile(fname)

print(f"Written {rank} face files")
print(f"File size each: {N_FIELDS * (sNx + 1) * (sNy + 1) * 8} bytes")
