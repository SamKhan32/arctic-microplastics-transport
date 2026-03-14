# extract_arctic_samples.py
# Run from shared/preprocessing/
# Extracts Arctic samples (lat >= 60N) to a clean CSV

import pandas as pd

XLSX = "original_data/marine_microplastics_samples.xlsx"
OUT  = "original_data/arctic_mp_samples.csv"
ARCTIC_LAT = 60.0

COLS = [
    "Record",
    "Sample Latitude",
    "Sample Longitude",
    "Sample Location",
    "Sample Year",
    "Ocean basin",
    "Total Pieces",
    "Total Pieces/L",
]

df = pd.read_excel(XLSX)

has_coords = df["Sample Latitude"].notna() & df["Sample Longitude"].notna()
arctic = df[has_coords & (df["Sample Latitude"] >= ARCTIC_LAT)][COLS].copy()

arctic.to_csv(OUT, index=True, index_label="original_index")

print(f"Wrote {len(arctic)} records to {OUT}")