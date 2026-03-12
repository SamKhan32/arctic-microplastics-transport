# NOTE: UNFINISHED
Currently, this baseline doesn't run - creating my own curvilinear proved to be a nightmare. The true first experiemnt is baseline_exp_llc270.
I do not guarantee that the code used here will work for anything.
# Arctic Circulation - Baseline Experiment

Baseline Arctic Ocean simulation using MITgcm with wind-only forcing on a curvilinear grid.

## Grid

| Parameter | Value |
|-----------|-------|
| Pole location | (-40E, 75N) - interior Greenland |
| Dimensions | 1530 x 520 cells |
| Resolution | ~10.7 km |
| Coverage | Full Arctic Ocean + key straits |
| Domain | Rotated lon -90 to 90, rotated lat 40 to 90 |

Curvilinear displaced-pole grid. Closed basin, no open boundaries.

## Physics

Wind-only forcing, barotropic, no tracers. Intended as a baseline to isolate the Arctic Ocean's response to climatological wind stress before adding thermodynamics or baroclinic structure.

## Input Code

All preprocessing scripts live in `input_code/`. Outputs go to `mitgcm_input/` and diagnostic plots to `analysis/`.

### `config.py`
Shared grid parameters and `make_grid()` function. All other scripts import from here. Defines the displaced-pole rotation and returns `LON`, `LAT` as `(NY, NX)` geographic coordinate arrays.

### `bathy.py`
Takes IBCAO v5 (500m resolution) as input. Sets all depths above sea level to zero, interpolates onto the curvilinear grid (MITgcm expects depths in meters, negative down), flood-fills to isolate the Arctic basin, and manually closes the Bering Strait. Writes `mitgcm_input/bathy.bin` and saves a diagnostic plot to `analysis/`.

### `metrics.py`
Computes MITgcm metric files from the `LON`/`LAT` arrays using finite differencing on the sphere. Writes the following big-endian float32 binaries to `mitgcm_input/`:

- `DXG.bin`, `DYG.bin` - grid face lengths [m]
- `DXC.bin`, `DYC.bin` - cell center spacings [m]
- `RAC.bin` - cell areas [m^2]
- `angleCosC.bin`, `angleSinC.bin` - grid orientation angle relative to true east

### `windStress.py`
Downloads ERA5 1991-2020 monthly mean surface stress (`avg_iews`, `avg_inss`) via the CDS API, averages into a 12-month climatology, interpolates onto the curvilinear grid, rotates from geographic (east/north) into grid-relative coordinates using the angle files, and masks to ocean points using the bathymetry. Writes `oceTauX_01.bin` through `oceTauX_12.bin` and `oceTauY_01.bin` through `oceTauY_12.bin` to `mitgcm_input/`.

Run with `--no-download` to skip the CDS request and reuse an existing NetCDF.

### `visualize_grid.py`
Plots the curvilinear grid on an Arctic polar stereographic projection. Diagnostic only.

### `visualize_wind_stress.py`
Plots climatological wind stress vectors overlaid on bathymetry for one or all months. Vectors are rotated back to geographic coordinates for display. Run with `--month N` to plot a single month, or without arguments for all 12. Saves PNGs to `analysis/`.

## Original Data

Raw source files in `original_data/`:

- `ibcao_v5_1_2025_400m_ice.tif` - IBCAO v5 bathymetry (input to `bathy.py`)
- `era5_stress_monthly_1991_2020.nc` - ERA5 surface stress download (input to `windStress.py`)

## Dependencies

```
numpy
scipy
matplotlib
cartopy
netCDF4
cdsapi
```
