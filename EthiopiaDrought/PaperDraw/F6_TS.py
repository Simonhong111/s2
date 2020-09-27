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


fig = plt.figure(figsize=(11, 6))
plt.xticks([])
plt.yticks([])
plt.axis('off')

ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("Bega")
ax1.plot(range(2007,2019),schrf,'-*',label="chiprs p")
ax1.plot(range(2007,2019),scmrf,'-*',label="cmip5 p")

ax1.legend()
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("Bega")
ax2.plot(range(2007,2019),schpvi,'-*',label="chiprs pvi")
ax2.plot(range(2007,2019),scmpvi,'-*',label="cmip5 pvi")
ax2.legend()
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("Krmet")
ax3.plot(range(2007,2019),Lchrf,'-*',label="chiprs p")
ax3.plot(range(2007,2019),Lcmrf,'-*',label="cmip5 p")
ax3.legend()
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("Krmet")
ax4.plot(range(2007,2019),Lchpvi,'-*',label="chiprs pvi")
ax4.plot(range(2007,2019),Lcmpvi,'-*',label="cmip5 pvi")
ax4.legend()
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.2)

plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_6.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()


