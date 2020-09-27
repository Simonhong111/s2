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

def ChirpsMask(CHIRPSDir,SeasonType="Short"):

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

def CMIP5TS(CMIP5Dir,CHIRPSDir,SeasonType="Short"):

    Mask = ChirpsMask(CHIRPSDir,SeasonType)

    MutYearAve = []

    for year in range(2007,2019):

        cm_file = os.path.join(CMIP5Dir,
                                   "cmip5_{}_{}.tif".format(SeasonType, str(year)))
        cm_data = gdal.Open(cm_file).ReadAsArray()

        MutYearAve.append(cm_data[Mask != -9999].mean())

    return MutYearAve

def ChirpsTS(CHIRPSDir,SeasonType="Short"):
    Mask = ChirpsMask(CHIRPSDir, SeasonType)

    MutYearAve = []

    for year in range(2007,2019):

        ch_file = os.path.join(CHIRPSDir,
                                   "chirps-v2.0.{}_{}.tif".format(SeasonType, str(year)))

        ch_data = gdal.Open(ch_file).ReadAsArray()

        MutYearAve.append(ch_data[Mask != -9999].mean())

    return MutYearAve

def CMIP5PVITS(CMPDir,CHIRPSDir,SeasonType="Short"):
    Mask = ChirpsMask(CHIRPSDir, SeasonType)

    MutYearAve = []

    for year in range(2007,2019):

        cmp_file = os.path.join(CMPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        cmp_data = gdal.Open(cmp_file).ReadAsArray()

        MutYearAve.append(cmp_data[Mask != -9999].mean())

    return MutYearAve
def ChirpsPVITS(CHPDir,CHIRPSDir,SeasonType="Short"):
    Mask = ChirpsMask(CHIRPSDir, SeasonType)
    MutYearAve = []

    for year in range(2007,2019):

        chp_file = os.path.join(CHPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        chp_data = gdal.Open(chp_file).ReadAsArray()

        MutYearAve.append(chp_data[Mask != -9999].mean())

    return MutYearAve



yy = "2006-2018"
mm = 'Short'

cm_rf_path = r"D:\Cornell\EthiopianDrought\CMIPMonth\Big"
cm_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CMIP"
ch_rf_path = r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth\Big"
ch_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CHIRPS"


schrf = ChirpsTS(ch_rf_path,'Short')
schpvi = ChirpsPVITS(ch_pvi_path,ch_rf_path,'Short')
scmrf = CMIP5TS(cm_rf_path,ch_rf_path,'Short')
scmpvi = CMIP5PVITS(cm_pvi_path,ch_rf_path,'Short')

Lchrf = ChirpsTS(ch_rf_path,"Long")
Lchpvi = ChirpsPVITS(ch_pvi_path,ch_rf_path,"Long")
Lcmrf = CMIP5TS(cm_rf_path,ch_rf_path,"Long")
Lcmpvi = CMIP5PVITS(cm_pvi_path,ch_rf_path,"Long")

fig = plt.figure(figsize=(10,6))

spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])




ax1.text(0.5, 0.95,'(a)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes)
ax1.plot(range(2007,2019),schrf,'-*',label="chiprs p")
ax1.plot(range(2007,2019),scmrf,'-*',label="cmip5 p")
ax1.set_ylabel("Belg")
ax1.legend()
# ax1.set_xlabel("P(mm)",fontdict={'fontname':'Times New Roman','fontsize':16})
ax2.text(0.5, 0.95,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)
ax2.plot(range(2007,2019),schpvi,'-*',label="chiprs pvi")
ax2.plot(range(2007,2019),scmpvi,'-*',label="cmip5 pvi")
ax2.legend()
# ax2.set_xlabel("PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
ax3.text(0.5, 0.95,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax3.transAxes)
ax3.plot(range(2007,2019),Lchrf,'-*',label="chiprs p")
ax3.plot(range(2007,2019),Lcmrf,'-*',label="cmip5 p")
ax3.legend()
ax3.set_ylabel("Kiremit",fontdict={'fontname':'Times New Roman','fontsize':16})
ax3.set_xlabel("P(mm)",fontdict={'fontname':'Times New Roman','fontsize':16})

ax4.text(0.5, 0.95,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax4.transAxes)
ax4.plot(range(2007,2019),Lchpvi,'-*',label="chiprs pvi")
ax4.plot(range(2007,2019),Lcmpvi,'-*',label="cmip5 pvi")
ax4.legend()
ax4.set_xlabel("PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.2)

plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_6.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()


