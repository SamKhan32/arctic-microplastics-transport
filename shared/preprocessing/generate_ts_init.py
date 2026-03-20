import numpy as np
import os

# Idealized Arctic T/S profile (5 levels)
# Level definitions match delR=50,100,250,600,4000m
T_profile = [-1.8, -1.5,  0.5,  1.5, -0.5]
S_profile = [31.0, 33.0, 34.5, 34.9, 34.94]

NR = 5
NY = 270
NX = 270

# Shape: (Nr, Ny, Nx) -- MITgcm expects level-by-level, row-major
theta = np.zeros((NR, NY, NX), dtype='>f8')
salt  = np.zeros((NR, NY, NX), dtype='>f8')

for k in range(NR):
    theta[k, :, :] = T_profile[k]
    salt[k, :, :]  = S_profile[k]

out_dir = os.path.join('inputs')

theta.tofile(os.path.join(out_dir, 'theta_init.bin'))
salt.tofile(os.path.join(out_dir, 'salt_init.bin'))

print('wrote inputs/theta_init.bin  shape=%s  dtype=>f8' % str(theta.shape))
print('wrote inputs/salt_init.bin   shape=%s  dtype=>f8' % str(salt.shape))
print('theta range: %.2f to %.2f' % (theta.min(), theta.max()))
print('salt  range: %.2f to %.2f' % (salt.min(),  salt.max()))