from osgeo import gdal,osr,ogr
import glob,os

import numpy as np




modis = gdal.Open(r"D:\Cornell\MOD13A3\MOD13A3.A2000032.h21v08.006.2015138123529.hdf")
subdateset = modis.GetSubDatasets()[1][0]
input_shape = r'D:\Cornell\Liuyanyan\WaterShed\watshed.shp'
print(subdateset)
# ds = gdal.Warp(
#                    r"D:\Cornell\MOD13A1.A2017353.h21v07.006.2018004225132clip3.tif",
#                     subdateset,
#                    format='GTiff',
#                    # cutlineDSName=input_shape,  # or any other file format
#                    # cutlineDSName=None,
#                    # cutlineWhere="FIELD = 'whatever'",
#                    #  cropToCutline=True,
#                    # optionally you can filter your cutline (shapefile) based on attribute values
#                    )  # select the no data value you like
#
#







# spat = osr.SpatialReference()
# spat.ImportFromEPSG(4326)
# print(spat.CloneGeogCS())
# raster = gdal.Open(r"C:\Users\zmhwh\Downloads\MOD13A1.A2017353.h21v07.006.2018004225132.tif")
# mprj = raster.GetProjection()
# print(raster.GetGeoTransform())
# dst = osr.SpatialReference()
# dst.ImportFromWkt(str(mprj))
# ct = osr.CoordinateTransformation(spat,dst)
# a = [37.77640348543618, 10.499059883146908]
# x,y = ct.TransformPoint(a[1],a[0])[0:2]
# print(x,y)
# dat = raster.ReadAsArray().astype(np.float)
# from matplotlib import  pyplot as plt
# g =raster.GetGeoTransform()
# XX =(x-g[0])/g[1]
# YY = (g[3]-y)/g[1]
# print(XX,YY)
# plt.imshow(dat)
# plt.scatter(XX,YY)
# plt.show()

# print(dat[int(YY)][int(XX)])
input_shape = r'D:\Cornell\Liuyanyan\WaterShed\watshed.shp'

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(input_shape, 0)
layer = dataSource.GetLayer()
for feature in layer:

    # print("geo",feature.GetGeometryRef())
    print(feature.GetGeometryRef().GetSpatialReference())
