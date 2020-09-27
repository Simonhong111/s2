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

Newdirname = r'D:\Cornell\NewSIF005Clip'
GOdirname = r'D:\Cornell\GOSIFClip'
Edirname = r'D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia'
GO2dirname = r'D:\Cornell\GOSIFV002Clip'
MON =['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


def SIF(dirname,month,product):
    Value = []
    Tm = []
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        if product == 'NewSIF':
            file = os.path.join(dirname,"SIF005_{}{}.nc.tif".format(str(dt.year),str(dt.month).zfill(2)))
        elif product == "GOSIF":
            file = os.path.join(dirname, "GOSIF_{}.M{}.tif".format(str(dt.year), str(dt.month).zfill(2)))
        elif product == 'EVI':
            file = os.path.join(dirname, "{}.{}.01.tif".format(str(dt.year), str(dt.month).zfill(2)))
        print(file)
        Tm.append(dt.year)
        sif = gdal.Open(file).ReadAsArray()
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip','MCD12C1.A{}001.*.tif'.format(str(dt.year))))[0]
        crop = gdal.Open(cropland).ReadAsArray()
        sif[crop != 12] = -9999
        mean = sif[(sif >-999) & (sif <32766)].mean()
        # print(dt, sif[sif>=-999].std())
        Value.append(mean)
    return np.array(Value),np.array(Tm)

def fit(y):
    return signal.detrend(y)



fig = plt.figure(figsize=(16, 20))
plt.title("Time Series of Anomaly over Cropland\n", fontsize=20)

for month in range(11,13):
    ax = fig.add_subplot(2, 1, month-11+1)
    NSIF,NTm = SIF(Newdirname,month,product='NewSIF')
    # GSIF,GTm = SIF(GOdirname,month,product='GOSIF')
    GSIF2,GTm2 = SIF(GO2dirname,month,product='GOSIF')
    Evi,ETm = SIF(Edirname,month,product='EVI')
    #
    NSIF = fit(NSIF)
    # GSIF = fit(GSIF)
    GSIF2 = fit(GSIF2)
    # Evi = fit(Evi)

    NSIF = (NSIF-np.mean(NSIF))/np.std(NSIF)
    # GSIF = (GSIF-np.mean(GSIF))/np.std(GSIF)
    GSIF2 = (GSIF2 - np.mean(GSIF2))/np.std(GSIF2)
    # NSIF = NSIF/np.std(NSIF)
    # GSIF = GSIF/np.std(GSIF)


    Evi = (Evi-np.mean(Evi))/np.std(Evi)
    ax.set_title(MON[month-1])


    ax.plot(NTm,NSIF,c='r',label='NewSIF')
    ax.scatter(NTm,NSIF,c='r',s=10)
    # plt.plot(GTm,GSIF,c='g',label='GOSIF')
    # plt.scatter(GTm,GSIF,c='g')
    ax.plot(GTm2,GSIF2,c='y',label='GOSIF2')
    ax.scatter(GTm2,GSIF2,c='y',s=10)
    ax.plot(ETm,Evi,c='b',label='EVI')
    ax.scatter(ETm,Evi,c='b',s=10)
    ax.plot(GTm2,[0]*len(list(GTm2)),'k')
    ax.legend()
plt.show()
