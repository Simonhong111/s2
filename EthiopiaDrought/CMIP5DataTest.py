from osgeo import gdal,osr,ogr

import numpy as np
from netCDF4 import Dataset
# Geo: 1°55'17.30"N,129°24'1.77"E
# MGRS: 52NEH4454012386
# Proj: Geographic Lat/Lon, WGS-84
#
# cmip5_20060101.tif
#   File: 68.5136,69.9628
#   Data: [0.875886]

raster = gdal.Open(r"D:\Cornell\EthiopianDrought\CMIP5Daily\cmip5_20060101.tif")
geo_t = raster.GetGeoTransform()
data = raster.ReadAsArray()
print(data[69,68])
lon = geo_t[0] + 68.5136 *geo_t[1]
lat = geo_t[3] + 69.9628*geo_t[5]

path =r"C:\Users\zmhwh\Downloads\pr_day_ACCESS1-0_historicalExt_r2i1p1_20060101-20201231.nc"
#
fin = Dataset(path, "r")

pr = fin.variables["pr"][:]*3600*24
lat_bnds = fin.variables["lat_bnds"][:]
lon_bnds = fin.variables["lon_bnds"][:]
lat_s = lat_bnds[:,0]
lat_e = lat_bnds[:,1]
lon_s = lon_bnds[:,0]
lon_e = lon_bnds[:,1]

mask_lat = (lat >= lat_s) & (lat < lat_e)
mask_lon = (lon >= lon_s) & (lon < lon_e)
mask_lat = np.where(mask_lat)
mask_lon = np.where(mask_lon)
print(mask_lat)
print(mask_lon)

print(pr[0][mask_lat[0],mask_lon[0]])


import xarray as xr
# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

dset = xr.open_dataset(path)
clim = dset['pr'].mean('time', keep_attrs=True)
print("**",clim)