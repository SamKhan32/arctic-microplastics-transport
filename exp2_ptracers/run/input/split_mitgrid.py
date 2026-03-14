import numpy as np

Nx, Ny = 270, 270
sNx, sNy = 27, 30
nPx, nPy = 10, 9
n_fields = 16

data = np.fromfile('arctic_cap_v2.mitgrid', dtype='>f8')
fields = data.reshape(n_fields, Ny, Nx)

rank = 0
for j in range(nPy):
    for i in range(nPx):
        rank += 1
        x0, x1 = i*sNx, (i+1)*sNx
        y0, y1 = j*sNy, (j+1)*sNy

        # Extract sNx+1 x sNy+1, wrapping at boundaries
        x1e = min(x0+sNx+1, Nx)
        y1e = min(y0+sNy+1, Ny)

        tile = np.zeros((n_fields, sNy+1, sNx+1), dtype='>f8')
        tile[:, :y1e-y0, :x1e-x0] = fields[:, y0:y1e, x0:x1e]

        fname = f'arctic_cap_v2.mitgrid.face{rank:03d}.bin'
        tile.astype('>f8').tofile(fname)

print(f"Written {rank} face files")
