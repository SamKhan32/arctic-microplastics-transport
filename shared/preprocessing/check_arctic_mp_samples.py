# check_arctic_mp_samples.py
# Run from shared/preprocessing/
# Reports how many samples fall north of 60N (Arctic threshold)

import pandas as pd

XLSX = "original_data/marine_microplastics_samples.xlsx"
LAT_COL = "Sample Latitude"
LON_COL = "Sample Longitude"
ARCTIC_LAT = 60.0

df = pd.read_excel(XLSX)

total = len(df)
has_coords = df[LAT_COL].notna() & df[LON_COL].notna()
arctic = df[has_coords & (df[LAT_COL] >= ARCTIC_LAT)]

print(f"Total records:        {total}")
print(f"Records with coords:  {has_coords.sum()}")
print(f"Arctic (lat >= 60N):  {len(arctic)}")
print()
print(arctic[[LAT_COL, LON_COL, "Sample Location", "Sample Year", "Ocean basin"]].to_string())