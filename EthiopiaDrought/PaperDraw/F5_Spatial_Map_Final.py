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
import matplotlib.gridspec as gridspec

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

fraction =0.069
fig = plt.figure(figsize=(10,6))
spec = gridspec.GridSpec(ncols=4, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[0, 2])
ax4 = fig.add_subplot(spec[0, 3])
ax5 = fig.add_subplot(spec[1, 0])
ax6 = fig.add_subplot(spec[1, 1])
ax7 = fig.add_subplot(spec[1, 2])
ax8 = fig.add_subplot(spec[1, 3])



ax1.text(0.95, 0.05,'(a)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax1.transAxes)


cmrf[cmrf == -9999] = np.nan
cax1 = ax1.imshow(cmrf, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=6.5)
cbar1 = plt.colorbar(cax1, ax=ax1, fraction=fraction, pad=0.04)
ax1.plot(x,y)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_ylabel("Belg",fontdict={'fontname':'Times New Roman','fontsize':16})
print("d",np.array(x).min(),np.array(x).max())

# chiprs p

ax2.text(0.95, 0.05,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax2.transAxes)

chrf[chrf == -9999] = np.nan
cax2 = ax2.imshow(chrf, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=6.5)
print("chrf vmax",chrf.max())
cbar2 = plt.colorbar(cax2, ax=ax2, fraction=fraction, pad=0.04)
ax2.plot(x,y)
ax2.set_xticks([])
ax2.set_yticks([])
# cmip5 pvi

ax3.text(0.95, 0.05,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax3.transAxes)

cmpvi[cmpvi == -9999] = np.nan
cax3 = ax3.imshow(cmpvi, cmap=plt.get_cmap("RdBu"), vmin=0.15, vmax=0.7)
cbar3 = plt.colorbar(cax3, ax=ax3, fraction=fraction, pad=0.04)
ax3.plot(x,y)
ax3.set_xticks([])
ax3.set_yticks([])
# chirps pvi

ax4.text(0.95, 0.05,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax4.transAxes)

chpvi[chpvi == -9999] = np.nan
cax4 = ax4.imshow(chpvi, cmap=plt.get_cmap("RdBu"), vmin=0.15, vmax=0.7)
cbar4 = plt.colorbar(cax4, ax=ax4, fraction=fraction, pad=0.04)
ax4.plot(x,y)
ax4.set_xticks([])
ax4.set_yticks([])

#********************
yy = "2006-2018"
mm = 'Long'
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

ax5.text(0.95, 0.05,'(e)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax5.transAxes)


cmrf[cmrf == -9999] = np.nan
cax1 = ax5.imshow(cmrf, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=10)
cbar1 = plt.colorbar(cax1, ax=ax5, fraction=fraction, pad=0.04)
ax5.plot(x,y)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_ylabel("Kiremit",fontdict={'fontname':'Times New Roman','fontsize':16})
ax5.set_xlabel("Cmip5 P",fontdict={'fontname':'Times New Roman','fontsize':16})
print("d",np.array(x).min(),np.array(x).max())

# chiprs p

ax6.text(0.95, 0.05,'(f)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax6.transAxes)

chrf[chrf == -9999] = np.nan
cax2 = ax6.imshow(chrf, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=10)
print("chrf vmax",chrf.max())
cbar2 = plt.colorbar(cax2, ax=ax6, fraction=fraction, pad=0.04)
ax6.plot(x,y)
ax6.set_xticks([])
ax6.set_yticks([])
ax6.set_xlabel("Chirps P",fontdict={'fontname':'Times New Roman','fontsize':16})
# cmip5 pvi

ax7.text(0.95, 0.05,'(g)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax7.transAxes)

cmpvi[cmpvi == -9999] = np.nan
cax3 = ax7.imshow(cmpvi, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=0.7)
cbar3 = plt.colorbar(cax3, ax=ax7, fraction=fraction, pad=0.04)
ax7.plot(x,y)
ax7.set_xticks([])
ax7.set_yticks([])
ax7.set_xlabel("Cmip5 PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
# chirps pvi

ax8.text(0.95, 0.05,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax8.transAxes)

chpvi[chpvi == -9999] = np.nan
cax4 = ax8.imshow(chpvi, cmap=plt.get_cmap("RdBu"), vmin=0.1, vmax=0.7)
cbar4 = plt.colorbar(cax4, ax=ax8, fraction=fraction, pad=0.04)
ax8.plot(x,y)
ax8.set_xticks([])
ax8.set_yticks([])
ax8.set_xlabel("Chirps PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()



