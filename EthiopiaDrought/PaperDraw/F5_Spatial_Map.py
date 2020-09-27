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


def CMIP5Comp(CMIP5Dir,SeasonType="Short"):

    YearNum = 13
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0
    for year in range(2007,2019):

        cm_file = os.path.join(CMIP5Dir,
                                   "cmip5_{}_{}.tif".format(SeasonType, str(year)))
        cm_data = gdal.Open(cm_file).ReadAsArray()

        Multidarr[:, :, band_id] = cm_data
        Mask[cm_data == -9999] = -9999
        band_id += 1

    SeIMG = Multidarr.mean(axis=2)
    SeIMG[Mask == -9999] = -9999
    return SeIMG

def ChirpsComp(CHIRPSDir,SeasonType="Short"):

    YearNum = 13
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0

    for year in range(2007,2019):

        ch_file = os.path.join(CHIRPSDir,
                                   "chirps-v2.0.{}_{}.tif".format(SeasonType, str(year)))

        ch_data = gdal.Open(ch_file).ReadAsArray()

        Multidarr[:, :, band_id] = ch_data
        Mask[ch_data == -9999] = -9999
        band_id += 1

    SeIMG = Multidarr.mean(axis=2)
    SeIMG[Mask == -9999] = -9999
    return SeIMG

def CMIP5PVIComp(CMPDir,SeasonType="Short"):

    YearNum = 13
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros((20, 13, YearNum))
    band_id = 0

    for year in range(2007,2019):

        cmp_file = os.path.join(CMPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        cmp_data = gdal.Open(cmp_file).ReadAsArray()

        Multidarr[:, :, band_id] = cmp_data
        Mask[cmp_data == -9999] = -9999
        band_id += 1

    SeIMG = Multidarr.mean(axis=2)
    SeIMG[Mask ==-9999] = -9999
    return SeIMG
def ChirpsPVIComp(CHPDir,SeasonType="Short"):

    YearNum = 13
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0

    for year in range(2007,2019):

        chp_file = os.path.join(CHPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        chp_data = gdal.Open(chp_file).ReadAsArray()

        Multidarr[:, :, band_id] = chp_data
        Mask[chp_data == -9999] = -9999
        band_id += 1

    SeIMG = Multidarr.mean(axis=2)
    SeIMG[Mask ==-9999] = -9999
    return SeIMG




ref_path = r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI\Big\long_pvi_2006.tif"
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



yy = "2006-2018"
mm = 'Short'
cm_rf_path = r"D:\Cornell\EthiopianDrought\CMIPMonth\Big"
cm_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CMIP"
ch_rf_path = r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth\Big"
ch_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CHIRPS"


chrf = ChirpsComp(ch_rf_path,mm)
chpvi = ChirpsPVIComp(ch_pvi_path,mm)
cmrf = CMIP5Comp(cm_rf_path,mm)
cmpvi = CMIP5PVIComp(cm_pvi_path,mm)

cmrf[chrf==-9999]=-9999
cmpvi[chpvi==-9999]=-9999



print(cmrf.max(),cmrf[chrf >-9999].min())
print(chrf.max(),chrf[chrf >-9999].min())
print(cmpvi.max(),cmpvi[chrf >-9999].min())
print(chpvi.max(),chpvi[chpvi>-9999].min())

chrfmax,chrfmin = chrf.max(),chrf[chrf >-9999].min()
chpvimax,chpvimin = chpvi.max(),chpvi[chpvi>-9999].min()
cmrfmax,cmrfmin = cmrf.max(),cmrf[chrf >-9999].min()
cmpvimax,cmpvimin = cmpvi.max(),cmpvi[chrf >-9999].min()


fig = plt.figure(figsize=(4, 5))
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.text(0.05, 0.95,'(e)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)

mask1 = np.where(cmrf > -9999)
cmrf[cmrf == -9999] = np.nan
cax1 = plt.imshow(cmrf, cmap=plt.get_cmap("RdBu"), vmin=0, vmax=10)
cbar1 = plt.colorbar(cax1, ax=plt.gca(), fraction=0.07, pad=0.04)
plt.plot(x,y)
print("d",np.array(x).min(),np.array(x).max())
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5a.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# chiprs p
fig = plt.figure(figsize=(4, 5))
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.text(0.05, 0.95,'(f)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
mask2 = np.where(chrf > -9999)
chrf[chrf == -9999] = np.nan
cax2 = plt.imshow(chrf, cmap=plt.get_cmap("RdBu"), vmin=0, vmax=10)
print("chrf vmax",chrf.max())
cbar2 = plt.colorbar(cax2, ax=plt.gca(), fraction=0.07, pad=0.04)
plt.plot(x,y)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5b.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# cmip5 pvi
fig = plt.figure(figsize=(4, 5))
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.text(0.05, 0.95,'(g)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
mask3 = np.where(cmpvi > -9999)
cmpvi[cmpvi == -9999] = np.nan

cax3 = plt.imshow(cmpvi, cmap=plt.get_cmap("RdBu"), vmin=0.15, vmax=0.7)
cbar3 = plt.colorbar(cax3, ax=plt.gca(), fraction=0.07, pad=0.04)
plt.plot(x,y)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5c.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# chirps pvi
fig = plt.figure(figsize=(4, 5))
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.text(0.05, 0.95,'(h)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
mask4 = np.where(chpvi > -9999)
chpvi[chpvi == -9999] = np.nan
cax4 = plt.imshow(chpvi, cmap=plt.get_cmap("RdBu"), vmin=0.15, vmax=0.7)
cbar4 = plt.colorbar(cax4, ax=plt.gca(), fraction=0.07, pad=0.04)
plt.plot(x,y)

plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5d.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()


