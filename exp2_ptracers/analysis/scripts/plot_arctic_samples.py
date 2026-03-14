# plot_arctic_samples.py
# Run from exp2_ptracers/
# Plots the 68 Arctic microplastic sample locations on a polar stereographic projection

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

CSV = "../shared/original_data/arctic_mp_samples.csv"

df = pd.read_csv(CSV)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.OCEAN, color='lightblue')
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.gridlines(draw_labels=False, linewidth=0.4, color='gray', linestyle='--')

ax.scatter(
    df["Sample Longitude"],
    df["Sample Latitude"],
    transform=ccrs.PlateCarree(),
    s=20,
    color='red',
    zorder=5,
    label=f"Samples (n={len(df)})"
)

ax.set_title("Arctic Microplastic Sample Locations")
plt.legend(loc='lower left')
plt.tight_layout()
plt.savefig("analysis/figures/prerun/arctic_samples.png", dpi=150)
print("Saved to analysis/figures/prerun/arctic_samples.png")