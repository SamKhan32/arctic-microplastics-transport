import netCDF4 as nc
ds = nc.Dataset(r'C:\Users\psam1\Desktop\ourProjects\arctic_circulation_spring\baseline_exp_llc270\original_data\ECCO-GRID_06.nc')
print(list(ds.variables.keys()))
ds.close()