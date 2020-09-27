from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import signal

def clipbyshp(input_raster,output_raster,input_shape, dstNodata=-9999):
    """
    :param input_raster: the raster data being processed later
    :param output_raster: the clipped datas' savepaths
    :param input_shape: the shape defining the extent
    :return: none
    """
    ds = gdal.Warp(output_raster,
                   input_raster,
                   format='GTiff',
                   cutlineDSName=input_shape,  # or any other file format
                   # cutlineDSName=None,
                   # cutlineWhere="FIELD = 'whatever'",
                   # optionally you can filter your cutline (shapefile) based on attribute values
                   cropToCutline=True,
                   dstNodata=dstNodata)  # select the no data value you like
    ds = None

def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset



def chirpsAnomMap(chirpsdirectory,yy,mm):


    chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}_{}.tif".format(mm,str(yy)))
    chirps_raster = gdal.Open(chirps_file).ReadAsArray()

    return chirps_raster

# anomalyMap = chirpsAnomMap(r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth",2009,'short')
# anomalyMap[anomalyMap==-9999] = np.nan
# plt.imshow(anomalyMap)
# plt.colorbar()
# plt.show()


def chirpsPVIAnomMap(chirpsdirectory,yy,mm):


    chirps_file = os.path.join(chirpsdirectory,"{}_pvi_{}.tif".format(mm,str(yy)))
    chirps_raster = gdal.Open(chirps_file).ReadAsArray()

    return chirps_raster



def cmip5AnomMap(chirpsdirectory,yy,mm):


    chirps_file = os.path.join(chirpsdirectory,"cmip5_{}_{}.tif".format(mm,str(yy)))
    chirps_raster = gdal.Open(chirps_file).ReadAsArray()

    return chirps_raster

# anomalyMap = chirpsAnomMap(r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth",2009,'short')
# anomalyMap[anomalyMap==-9999] = np.nan
# plt.imshow(anomalyMap)
# plt.colorbar()
# plt.show()


def cmip5PVIAnomMap(chirpsdirectory,yy,mm):


    chirps_file = os.path.join(chirpsdirectory,"{}_pvi_{}.tif".format(mm,str(yy)))
    chirps_raster = gdal.Open(chirps_file).ReadAsArray()


    return chirps_raster






# cm_rf = gdal.Open(r"D:\Cornell\EthiopianDrought\CMIPMonth\cmip5_long_2006.tif")
# cm_pvi = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI\long_pvi_2006.tif")
# ch_rf = gdal.Open(r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth\agg_chirps-v2.0.long_2006.tif")
# ch_pvi = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVIDaily\long_pvi_2006.tif")
#
#
#
#
# print(cm_rf.GetGeoTransform())
# print(cm_pvi.GetGeoTransform())
# print(ch_rf.GetGeoTransform())
# print(ch_pvi.GetGeoTransform())


# 参考


ref_path = r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI\long_pvi_2006.tif"
ref_raster = gdal.Open(ref_path)
geo_t = ref_raster.GetGeoTransform()

# 计算矢量边界

daShapefile = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()
geo = str(geo).split("((")[1].split("))")[0].split(",")
x = []
y = []
for term in geo:
    x.append(float(term.split(" ")[0]))
    y.append(float(term.split(" ")[1]))

x = np.array(x)
y = np.array(y)
x = (x - geo_t[0]) / geo_t[1]
y = (y - geo_t[3]) / geo_t[5]

# plt.imshow(ref_raster.ReadAsArray())
# # plt.colorbar()
# plt.plot(x,y)
# plt.show()


yy = str(2015)
mm = 'long'
cm_rf_path = r"D:\Cornell\EthiopianDrought\CMIPMonth"
cm_pvi_path = r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI"
ch_rf_path = r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth"
ch_pvi_path = r"D:\Cornell\EthiopianDrought\AData\PVIDaily"


chrf = chirpsAnomMap(ch_rf_path,yy,mm)
chpvi = chirpsPVIAnomMap(ch_pvi_path,yy,mm)
cmrf = cmip5AnomMap(cm_rf_path,yy,mm)
cmpvi = cmip5PVIAnomMap(cm_pvi_path,yy,mm)

print(cmrf.max(),cmrf.min())
print(chrf.max(),chrf[chrf >-9999].min())
print(cmpvi.max(),cmpvi.min())
print(chpvi.max(),chpvi[chpvi>-9999].min())

chrfmax,chrfmin = chrf.max(),chrf[chrf >-9999].min()
chpvimax,chpvimin = chpvi.max(),chpvi[chpvi>-9999].min()
cmrfmax,cmrfmin = cmrf.max(),cmrf.min()
cmpvimax,cmpvimin = cmpvi.max(),cmpvi.min()


fig = plt.figure(figsize=(10, 15))
plt.title("{} {} rains Original Value Map\n ".format(yy,mm) + '\n', fontsize=16)
plt.xticks([])
plt.yticks([])

ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("{} Rains CMIP5 Rainfall Original Value Map".format(mm))
mask1 = np.where(cmrf > -9999)
cmrf[cmrf == -9999] = np.nan
cax1 = ax1.imshow(cmrf, cmap=plt.get_cmap("RdBu"), vmin=0, vmax=10)
cbar1 = plt.colorbar(cax1, ax=ax1, fraction=0.036, pad=0.04)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title("{} rains average (mm)".format(mm))
ax1.plot(x,y)

ax2 = fig.add_subplot(2, 2,2)
ax2.set_title("{} Rains CHIRPS Rainfall Original Value Map".format(mm))
mask2 = np.where(chrf > -9999)
chrf[chrf == -9999] = np.nan
cax2 = ax2.imshow(chrf, cmap=plt.get_cmap("RdBu"), vmin=0, vmax=10)
print("chrf vmax",chrf.max())
cbar2 = plt.colorbar(cax2, ax=ax2, fraction=0.036, pad=0.04)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_ylabel("{} rains average (mm)".format(mm))
ax2.plot(x,y)

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("{} Rains CMIP5 PVI Original Value Map".format(mm))
mask3 = np.where(cmpvi > -9999)
cmpvi[cmpvi == -9999] = np.nan

cax3 = ax3.imshow(cmpvi, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=1)
cbar3 = plt.colorbar(cax3, ax=ax3, fraction=0.036, pad=0.04,size=0.2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.plot(x,y)

ax4 = fig.add_subplot(2, 2,4)
ax4.set_title("{} Rains CHIRPS PVI Original Value Map".format(mm))
mask4 = np.where(chpvi > -9999)
chpvi[chpvi == -9999] = np.nan
cax4 = ax4.imshow(chpvi, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=1)
cbar4 = plt.colorbar(cax4, ax=ax4, fraction=0.036, pad=0.04)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.plot(x,y)
# fig.tight_layout()  # 调整整体空白
plt.show()


