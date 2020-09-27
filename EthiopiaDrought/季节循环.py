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

def EviM(evidir,cropdir,mm):

    croppath = os.path.join(cropdir,"agg_clip.tif")
    crop = gdal.Open(croppath).ReadAsArray()
    Average = []
    Num = []
    for yy in range(2003,2019):
        evipath = glob.glob(os.path.join(evidir, "{}.{}.01.tif".format(str(yy), str(mm).zfill(2))))[0]
        evi = gdal.Open(evipath).ReadAsArray()
        mask = (crop == 255) & (evi > -3000)
        evimask = evi[mask]
        sum = evimask.sum()*0.0001
        Average.append(sum)
        Num.append(evimask.shape[0])
        del evi
        del mask
        del evimask

    return np.array(Average).sum()/np.array(Num).sum()

def RFM(RFdir,cropdir,mm):
    croppath = os.path.join(cropdir, "agg_clip.tif")
    crop = gdal.Open(croppath).ReadAsArray()
    Average = []
    Num = []
    for yy in range(2003,2019):
        RFpath = glob.glob(os.path.join(RFdir,"chirps-v2.0.{}.{}.tif".format(str(yy),str(mm).zfill(2))))[0]
        RF = gdal.Open(RFpath).ReadAsArray()
        mask = (crop == 255) & (RF > -9999)
        RFmask = RF[mask]
        sum = RFmask.sum()
        Average.append(sum)
        Num.append(RFmask.shape[0])
        del RF
        del mask
        del RFmask


    return np.array(Average).sum()/np.array(Num).sum()

def NewSIF(NSIFdir,cropdir,mm):
    croppath = os.path.join(cropdir, "agg_clip.tif")
    crop = gdal.Open(croppath).ReadAsArray()
    Average = []
    Num = []
    for yy in range(2003, 2019):
        NSIFpath = glob.glob(os.path.join(NSIFdir, "SIF005_{}{}.nc.tif".format(str(yy), str(mm).zfill(2))))[0]
        NSIF = gdal.Open(NSIFpath).ReadAsArray()
        mask = (crop == 255) & (NSIF > -9999)
        NSIFmask = NSIF[mask]
        sum = NSIFmask.sum()
        Average.append(sum)
        Num.append(NSIFmask.shape[0])
        del NSIF
        del mask
        del NSIFmask
    return np.array(Average).sum() / np.array(Num).sum()


def GOSIF(GOSIFdir,cropdir,mm):
    croppath = os.path.join(cropdir, "agg_clip.tif")
    crop = gdal.Open(croppath).ReadAsArray()
    Average = []
    Num = []
    for yy in range(2003, 2019):
        GOSIFpath = glob.glob(os.path.join(GOSIFdir, "GOSIF_{}.M{}.tif".format(str(yy), str(mm).zfill(2))))[0]
        GOSIF = gdal.Open(GOSIFpath).ReadAsArray()
        mask = (crop == 255) & (GOSIF < 32766)
        GOSIFmask = GOSIF[mask]
        sum = GOSIFmask.sum()*0.0001
        Average.append(sum)
        Num.append(GOSIFmask.shape[0])
        del GOSIF
        del mask
        del GOSIFmask
    return np.array(Average).sum() / np.array(Num).sum()




evidir = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
RFdir =r"D:\Cornell\EthiopianDrought\Chirps2"
NSIFdir =r"D:\Cornell\NewSIF005Clip"
GOSIFdir =r"D:\Cornell\GOSIFV002Clip"
cropdir = r"D:\Cornell\EthiopianDrought\CropType2015"
EMean = []
RMean = []
NMean = []
GMean = []

for mm in range(1,13):
    EMean.append(EviM(evidir,cropdir,mm))
    RMean.append(RFM(RFdir,cropdir,mm))
    NMean.append(NewSIF(NSIFdir,cropdir,mm))
    GMean.append(GOSIF(GOSIFdir,cropdir,mm))

Month = range(1,13)
fig = plt.figure(figsize=(14, 12))
plt.title("mean seasonal cycle without detrend"+ '\n', fontsize=14)
plt.xticks([])
plt.yticks([])
ax = fig.add_subplot(2,2,1)
ax.set_title("evi")
ax.plot(Month,EMean,label="evi")
ax.scatter(Month,EMean)
ax.legend()

ax2 = fig.add_subplot(2,2,2)
ax2.set_title("rainfall")
ax2.set_ylabel("mm", fontsize=20)
ax2.plot(Month,RMean,label='rainfall')
ax2.scatter(Month,RMean)
ax2.legend()

ax3 = fig.add_subplot(2,2,3)
ax3.set_title("newsif")
ax3.set_ylabel(r'W$m^{-2} \mu^{-1} sr^{-1}$', fontsize=20)
ax3.plot(Month,NMean,label='newsif')
ax3.scatter(Month,NMean)
ax3.legend()
ax4 = fig.add_subplot(2,2,4)
ax4.set_title("gosif")
ax4.set_ylabel(r'W$m^{-2} \mu^{-1} sr^{-1}$', fontsize=20)
ax4.plot(Month,GMean,label='gosif')
ax4.scatter(Month,GMean)
ax4.legend()


# plt.legend()
plt.show()


