import os
import numpy
from GridGen import *
from osgeo import gdal,ogr,osr
shp_path = r"D:\Sen2Projecton\HBSHP\Hubei.shp"
# grid_path = r""
tile_path = r"D:\Sen2Projecton\HBS2TILE\HBS2.shp"

#
# newGridGen = GridDefination()
# mGrid = newGridGen.validGrid(shp_path,tile_path," ")

# import  time
# start = time.time()
#
# gridtilepairs = np.load(r"D:\Sen2Projecton\np.npz")
#
# end = time.time()
#
# print(end-start)
#
# print(gridtilepairs["gridId"])
# print(gridtilepairs["Tile"])
# print(gridtilepairs["EtRing"])






# maskband = gdal.Open(r"D:\Sen2Projecton\HBCrop\fromglc10v01_28_108.tif")
# help(maskband)
# print(maskband.GetProjection())
# # print(maskband.GetMetadata())
# # help(maskband)
# print(maskband.GetGeoTransform())
# GT = maskband.GetGeoTransform()
# maskraster = maskband.ReadAsArray().astype(np.float)
# cols = maskband.RasterXSize
# rows = maskband.RasterYSize
# Xgeo = GT[0] + cols*GT[1] + rows*GT[2]
# Ygeo = GT[3] + cols*GT[4] + rows*GT[5]
# print(Xgeo,Ygeo)
# lond = int(Xgeo)
# lonm = int((Xgeo-lond)*60)
# lons = (Xgeo - lond - lonm/60.0)*3600
#
# land = int(Ygeo)
# lanm = int((Ygeo-land)*60)
# lans = (Ygeo - land - lanm/60.0)*3600
#
# print(lond,lonm,lons)
# print(land,lanm,lans)








from FileteredGrid import *

# writeLandCoverNpz(r"D:\Sen2Projecton\HBCrop",r"D:\Sen2Projecton\landcover.npz")

valid_grid_path = r"D:\Sen2Projecton\np.npz"
landcover_path = r"D:\Sen2Projecton\HBCrop\HBCrop.tif"
savefilter_path = r"D:\Sen2Projecton\FilterGrid.npz"



# valGridId, valGridTile,valGridEtRing = writeFilterGrid(valid_grid_path,landcover_path,savefilter_path)
#
#
# filter(landcover_path,valGridId, valGridTile,valGridEtRing, r"D:\Sen2Projecton\Filter.npz")
filternpz_path = r"D:\Sen2Projecton\Filter.npz"
save_shp_path = r"D:\Sen2Projecton\HBSHP\FilterGrid.shp"
# saveFilterGrid2Shp(filternpz_path,save_shp_path)
sen2dir = r"J:\2018S2docker"

from Sen2Aggragation import *

# readSen2(filternpz_path,sen2dir,r"D:\Sen2Projecton\Sen2Mean.npz")

#
# path = r"J:\2018S2docker\T49RGN\S2A_MSIL2A_20180513T030541_N0206_R075_T49RGN_20180513T064347.SAFE"
#
# mSenL2A = SenL2AReader(path)
# cloudPro = mSenL2A.getQiData("MSK_CLDPRB")
# print(cloudPro.GetGeoTransform())
#
#
#
# prosrs = osr.SpatialReference()
# prosrs.ImportFromWkt(cloudPro.GetProjection())
# geosrs = prosrs.CloneGeogCS()
# # print(prosrs)
# # print(geosrs)
# ct = osr.CoordinateTransformation(geosrs, prosrs)
#
# coords = ct.TransformPoint(29.81426153422464, 113.06917639199416)
# print(coords)
#
#
# print(geo2lonlat(cloudPro,699960.0,3300000.0))
#
# print(geo2imagexy(cloudPro,699960.0,3300000.0))


# a = np.array([[1,2],[3,4]])
#
# w = a<1
# w2 = (a>4) | (a <0)
#
# w3 = np.where(w2 == True)
# if not np.empty(w3[0]):
#     print("**")
# print(w3[0])
# # print(w3)
#
# print(a[w3])


import sys
from optparse import OptionParser
import h5py
import time
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
import time
import glob
from numba import jit


start = '2018-05-01'
stop = '2018-05-31'



from SIFAggaragation import *

data = gdal.Open(r"C:\Users\Administrator\Desktop\L2A_T50SKA_20170724T025551_B04_20m.jp2")

help(data)
print(data.GetGeoTransform())



trans = data.GetGeoTransform()
col = data.RasterXSize
row = data.RasterYSize
print(col,row)
px = trans[0] + col * trans[1] + row * trans[2]
py = trans[3] + col * trans[4] + row * trans[5]


print(geo2lonlat(data,trans[0],trans[3]))
print(geo2lonlat(data,px,py))