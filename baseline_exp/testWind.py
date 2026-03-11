import numpy as np

NX, NY = 1530, 520

bathy = np.fromfile('mitgcm_input/bathy.bin', dtype='>f4').reshape(NY, NX)
taux  = np.fromfile('mitgcm_input/oceTauX.bin', dtype='>f4').reshape(12, NY, NX)
tauy  = np.fromfile('mitgcm_input/oceTauY.bin', dtype='>f4').reshape(12, NY, NX)

print(f"bathy  min={bathy.min():.1f}  max={bathy.max():.1f}")
print(f"taux   min={taux.min():.4f}  max={taux.max():.4f}")
print(f"tauy   min={tauy.min():.4f}  max={tauy.max():.4f}")

print(f"bathy NaN/Inf: {np.isnan(bathy).sum()} / {np.isinf(bathy).sum()}")
print(f"taux  NaN/Inf: {np.isnan(taux).sum()} / {np.isinf(taux).sum()}")
print(f"tauy  NaN/Inf: {np.isnan(tauy).sum()} / {np.isinf(tauy).sum()}")

ocean_mask = bathy < 0
mag = np.sqrt(taux[0]**2 + tauy[0]**2)

print(f"\nJanuary stress over ocean:")
print(f"  max magnitude = {mag[ocean_mask].max():.4f} N/m^2")
print(f"  mean magnitude = {mag[ocean_mask].mean():.4f} N/m^2")
print(f"  stress on land cells (should be 0): {mag[~ocean_mask].max():.6f}")
import numpy as np
bathy = np.fromfile('mitgcm_input/bathy.bin', dtype='>f4').reshape(NY, NX)  # or f8 if float64
wet = (bathy < 0).astype(int)

# Look for isolated wet cells (all 4 neighbors are land)
from scipy.ndimage import label
labeled, nfeatures = label(wet)
print(f"Number of disconnected wet regions: {nfeatures}")
# Should be 1 for a closed Arctic basin. If >1, you have isolated pockets.