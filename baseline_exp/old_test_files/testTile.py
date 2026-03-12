import numpy as np

NX, NY = 1530, 520
fields = ['xC','yC','dxF','dyF','rA','xG','yG','dxV','dyU','rAz',
          'dxC','dyC','rAw','rAs','dxG','dyG','angleCosC','angleSinC']

data = np.fromfile('mitgcm_input/tile001.mitgrid', dtype='>f8')
print(f"Total elements: {len(data)}  expected: {18*NX*NY}")

for i, name in enumerate(fields):
    arr = data[i*NX*NY:(i+1)*NX*NY].reshape(NY, NX)
    print(f"{name:12s}  min={arr.min():.4g}  max={arr.max():.4g}  zeros={(arr==0).sum()}  neg={(arr<0).sum()}")
arr_dxF = data[2*NX*NY:3*NX*NY].reshape(NY, NX)
arr_dxC = data[10*NX*NY:11*NX*NY].reshape(NY, NX)

print("dxF < 1m locations:", np.argwhere(arr_dxF < 1))
print("dxC == 0 locations:", np.argwhere(arr_dxC == 0)[:10])
arr_dxF = data[2*NX*NY:3*NX*NY].reshape(NY, NX)
print("dxF > 1e6 locations:", np.argwhere(arr_dxF > 1e6)[:10])
import matplotlib.pyplot as plt
import numpy as np

NX, NY = 1530, 520
data = np.fromfile('mitgcm_input/tile001.mitgrid', dtype='>f8')
dxC = data[10*NX*NY:11*NX*NY].reshape(NY, NX)

plt.figure(figsize=(15,5))
plt.pcolormesh(dxC/1000, vmin=0, vmax=15)
plt.colorbar(label='dxC (km)')
plt.title('dxC')
plt.savefig('dxC_check.png', dpi=100)
print(f"dxC < 5km: {(dxC < 5000).sum()} cells")
print(f"dxC < 1km: {(dxC < 1000).sum()} cells")
