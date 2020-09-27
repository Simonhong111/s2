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

def EviM(evidir,cropdir,yy,mm):
    evipath = glob.glob(os.path.join(evidir,"{}.{}.01.tif".format(str(yy),str(mm).zfill(2))))[0]
    croppath = glob.glob(os.path.join(cropdir,"MCD12C1.A{}001.006.*.hdf.tif").format(str(yy)))[0]
    print(croppath)
    evi = gdal.Open(evipath).ReadAsArray()
    crop = gdal.Open(croppath).ReadAsArray()
    mask = (crop == 12) & (evi > -3000)
    evimask = evi[mask]
    mean = evimask.mean()*0.0001
    return mean

def RFM(RFdir,cropdir,yy,mm):
    RFpath = glob.glob(os.path.join(RFdir,"chirps-v2.0.{}.{}.tif".format(str(yy),str(mm).zfill(2))))[0]
    croppath = glob.glob(os.path.join(cropdir,"MCD12C1.A{}001.006.*.hdf.tif").format(str(yy)))[0]
    RF = gdal.Open(RFpath).ReadAsArray()
    crop = gdal.Open(croppath).ReadAsArray()
    mask = (crop == 12) & (RF > -9999)
    RFmask = RF[mask]
    mean = RFmask.mean()
    return mean

def NewSIF(NSIFdir,cropdir,yy,mm):
    NSIFpath = glob.glob(os.path.join(NSIFdir,"SIF005_{}{}.nc.tif".format(str(yy),str(mm).zfill(2))))[0]
    croppath = glob.glob(os.path.join(cropdir,"MCD12C1.A{}001.006.*.hdf.tif").format(str(yy)))[0]
    NSIF = gdal.Open(NSIFpath).ReadAsArray()
    crop = gdal.Open(croppath).ReadAsArray()
    mask = (crop == 12) & (NSIF > -9999)
    NSIFmask = NSIF[mask]
    mean = NSIFmask.mean()
    return mean

def GOSIF(GOSIFdir,cropdir,yy,mm):
    GOSIFpath = glob.glob(os.path.join(GOSIFdir,"GOSIF_{}.M{}.tif".format(str(yy),str(mm).zfill(2))))[0]
    croppath = glob.glob(os.path.join(cropdir,"MCD12C1.A{}001.006.*.hdf.tif").format(str(yy)))[0]
    GOSIF = gdal.Open(GOSIFpath).ReadAsArray()
    crop = gdal.Open(croppath).ReadAsArray()
    mask = (crop == 12) & (GOSIF < 32766)
    GOSIFmask = GOSIF[mask]
    mean = GOSIFmask.mean()*0.0001
    return mean


evidir = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
RFdir =r"D:\Cornell\EthiopianDrought\Chirps2"
NSIFdir =r"D:\Cornell\NewSIF005Clip"
GOSIFdir =r"D:\Cornell\GOSIFV002Clip"
cropdir = r"D:\Cornell\MCD12C1V006Clip"
EMean = []
RMean = []
NMean = []
GMean = []
yy = 2017
for mm in range(1,13):
    EMean.append(EviM(evidir,cropdir,yy,mm))
    RMean.append(RFM(RFdir,cropdir,yy,mm))
    NMean.append(NewSIF(NSIFdir,cropdir,yy,mm))
    GMean.append(GOSIF(GOSIFdir,cropdir,yy,mm))

Month = range(1,13)
fig = plt.figure(figsize=(14, 12))
plt.title("{} mean seasonal cycle without detrend".format(yy)+ '\n', fontsize=14)
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


