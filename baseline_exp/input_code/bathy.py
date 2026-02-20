import numpy as np
import rasterio
from pathlib import Path

# anchor dir

EXP_DIR = Path(__file__).resolve().parents[1]
TIF_FILE_PATH = EXP_DIR._str + "/data/ibcao_v5_1_2025_400m_ice.tif"

def tif_to_mitgcm_binary(tif_path, out_path, dtype='>f4'):
    """
    Convert a GeoTIFF to a flat binary file in Fortran (column-major) order
    suitable for MitGCM.
    
    dtype: MitGCM expects big-endian 32-bit float by default ('>f4')
           Use '>f8' if you're running in double precision
    """
    with rasterio.open(tif_path) as src:
        data = src.read(1)  # Read first band -> shape (rows, cols) = (ny, nx)
        
        print(f"Shape (ny, nx): {data.shape}")
        print(f"CRS: {src.crs}")
        print(f"Resolution: {src.res}")
        print(f"NoData value: {src.nodata}")
        
        # Replace nodata with 0.0 (MitGCM convention for land)
        if src.nodata is not None:
            data = np.where(data == src.nodata, 0.0, data)
    
    # Cast and byte-swap to big-endian float32
    data = data.astype(dtype)
    
    # Write in Fortran order: transpose first, then write F-order
    # MitGCM reads (nx, ny) with ny varying fastest in Fortran — 
    # so we need to transpose to (nx, ny) then write
    data_T = data.T  # shape (nx, ny)
    
    # Writing with order='F' on the transposed array gives correct memory layout
    data_T.flatten(order='F').tofile(out_path)
    
    print(f"Written to {out_path}: {data_T.shape} as {dtype}")

if __name__ == "__main__":
    out_file = EXP_DIR._str + "/data/bathy.bin"
    tif_to_mitgcm_binary(TIF_FILE_PATH, out_file)
