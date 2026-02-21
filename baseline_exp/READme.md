For my baseline arctic expeirment's input and output codes.
This is a: 
curvilinear grid
closed basin
~10km x 10km resolution (interpolated from IBACO's 500m x 500m grid)
only wind forcing

Here's a descrption of what is suppose to happen to all the input files
GRID
Pole:        (-40°E, 75°N) — interior Greenland
Dimensions:  1530 × 520 cells
Resolution:  ~10.7 km
Coverage:    Full Arctic Ocean + key straits
Domain:      Rotated lon -90→90°, rotated lat 40→90°
BATHYMETRY
bathy.py is suppose to take in the IBACO data, make all depths > 0, equal to 0, and then make a plot after printing out some data.
It will then interpolate the current resolution to the desired resolution (the mitgcm expects units of meters.)

