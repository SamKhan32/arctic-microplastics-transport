import numpy as np
import netCDF4 as nc

GRID_PATH = 'original_data/ECCO-GRID_06.nc'
IN_DIR    = 'mitgcm_input'
OUT_DIR   = 'mitgcm_input'

# Load rotation angles
g  = nc.Dataset(GRID_PATH)
CS = np.array(g['CS'][:].squeeze())  # cos(angle), (270, 270)
SN = np.array(g['SN'][:].squeeze())  # sin(angle), (270, 270)
g.close()

# Read geographic-frame stress (output of generate_wind_stress.py)
tx_geo = np.fromfile(f'{IN_DIR}/taux_jan.bin', dtype='>f8').reshape(270, 270)
ty_geo = np.fromfile(f'{IN_DIR}/tauy_jan.bin', dtype='>f8').reshape(270, 270)

# Rotate to grid-relative frame
tx_grid =  tx_geo * CS + ty_geo * SN
ty_grid = -tx_geo * SN + ty_geo * CS

tx_grid.astype('>f8').tofile(f'{OUT_DIR}/taux_jan_rot.bin')
ty_grid.astype('>f8').tofile(f'{OUT_DIR}/tauy_jan_rot.bin')

print(f"taux_rot range: {tx_grid.min():.4f} to {tx_grid.max():.4f}")
print(f"tauy_rot range: {ty_grid.min():.4f} to {ty_grid.max():.4f}")
print("Wrote taux_jan_rot.bin, tauy_jan_rot.bin")
