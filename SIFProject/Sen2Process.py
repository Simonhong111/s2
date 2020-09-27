from osgeo import gdal,osr,ogr
import sys,os
from Sentinel2Reader import SenL2AReader
import numpy as np
from GridGen import *
from Sen2Clip import *
from geo2mapxy import *
def dms2dec(degree, minute, second, sphere="N"):
    m = minute / 60.0
    s = second / 3600.0
    return degree + m + s

path = r'D:\S2A_MSIL2A_20190417T033541_N0211_R061_T48SWC_20190417T073942.SAFE'

daShapefile = r"D:\Sentinel-2Data\WHSHAPE\whan.shp"

ppath = r"D:\S2A_MSIL2A_20190417T033541_N0211_R061_T48SWC_20190417T073942.SAFE\GRANULE\L2A_T48SWC_A019935_20190417T033835\IMG_DATA\R20m\s2clsscom.tif"


print(dms2dec(105,10,11.88),dms2dec(34,12,2.59))
print(dms2dec(105,12,30.64),dms2dec(34,10,10.17))
scope =(105.16996666666667, 34.200719444444445,105.20851111111111, 34.169491666666666)
# GetRoi(ppath,r"D:\S2A_MSIL2A_20190417T033541_N0211_R061_T48SWC_20190417T073942.SAFE\GRANULE\L2A_T48SWC_A019935_20190417T033835\IMG_DATA\R20m\test.tif",scope)


from GridGen import *
from Aggaragation import  *

griddef = GridDefination()

# gridarr = griddef.getGrid(daShapefile)


sumpath = r"D:\Sen2data"

# mAggr = Sen2Aggragation(years=["2017"],directory=sumpath)
# mAggr.findSubImg("2017")

from datetime import *
import time
from dateutil import rrule
start = datetime.strptime("2019-05-01","%Y-%m-%d").date()
end = datetime.strptime("2019-05-31","%Y-%m-%d").date()

lastDay = datetime(2014,1,1)
last = datetime(2014,1,2)
delt = last - lastDay
print(delt)

# end = datetime.strptime(lastDay,"%Y-%m-%d").date()

agg = Sen2Aggragation
# agg.findSubImg("2017")
from datetime import datetime


# x = np.arange(1,100,1)
# y = np.sin(x)
# kwargs = {"arg1": x, "arg2": x, "arg3": y}
# # >>> test_args_kwargs(**kwargs)
# np.savez(r"D:\Sen2data\np.npz",**kwargs)
#
# data = np.load(r"D:\Sen2data\np.npz")
# print(data['arg1'])
#



shp_path = r"D:\Sen2Projecton\HBSHP\Hubei.shp"
outShp = r"D:\Sen2Projecton\HBSHP\HubeiGrid.shp"
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(shp_path, 0)
# print(dataSource.GetTransform())
layer = dataSource.GetLayer()
spatialRef = layer.GetSpatialRef()
print("spatialref",spatialRef.GetAuthorityCode("GEOGCS"))

# help(spatialRef)



mGrid = GridDefination()

mGrid.saveGrid(outShp,shp_path)

#
# import matplotlib.pyplot as plt
# fig = plt.figure()
# for m in mGrid:
#     lon = []
#     lon.append(float(m['Extent'][0]))
#     lon.append(float(m['Extent'][2]))
#     lon.append(float(m['Extent'][2]))
#     lon.append(float(m['Extent'][0]))
#     lon.append(float(m['Extent'][0]))
#     lat = []
#     lat.append(float(m['Extent'][1]))
#     lat.append(float(m['Extent'][1]))
#     lat.append(float(m['Extent'][3]))
#     lat.append(float(m['Extent'][3]))
#     lat.append(float(m['Extent'][1]))
#     plt.plot(lon,lat)
# plt.show()


