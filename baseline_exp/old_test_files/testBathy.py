import numpy as np
bathy = np.fromfile('mitgcm_input/bathy.bin', dtype='>f4').reshape(NY, NX)  # or f8 if float64
wet = (bathy < 0).astype(int)

# Look for isolated wet cells (all 4 neighbors are land)
from scipy.ndimage import label
labeled, nfeatures = label(wet)
print(f"Number of disconnected wet regions: {nfeatures}")
# Should be 1 for a closed Arctic basin. If >1, you have isolated pockets.