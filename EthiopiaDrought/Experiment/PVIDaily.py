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

def write_Img(data, path, im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    rainfall = gdal.Open(r"D:\Cornell\EthiopianDrought\ChirpsDaily2\chirps-v2.0.2003.01.01.tif")
    geotrans = rainfall.GetGeoTransform()
    proj = rainfall.GetProjection()
    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))


    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

def PVI(chirpsdirectory,year):
    start = datetime.strptime("-".join([str(year), str(2).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("-".join([str(year), str(5).zfill(2), "31"]), "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0.{}.{}.{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1
    print("bandNum",bandNum)
    mask = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):
        print(dt)
        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0.{}.{}.{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))

        chirps = gdal.Open(chirps_file).ReadAsArray()
        print(chirps_file,chirps[157,215])
        multidarr[:,:,band_id] = chirps
        mask[chirps == -9999] += 1
        del chirps
        band_id += 1
        axisTime.append([dt.year,dt.month])

    h, w, channels = multidarr.shape
    CC = np.zeros_like(multidarr,dtype=np.float)
    PP = np.zeros_like(multidarr,dtype=np.float)
    ppm = np.mean(multidarr, axis=2)
    for mc in range(1,channels+1):

        CC[:,:,mc-1] = np.sum(multidarr[:,:,0:mc],axis=2)


    EE = np.zeros_like(multidarr,dtype=np.float)
    for i in range(channels):
        EE[:,:,i] = ppm*(i+1)
    RR = CC/EE
    del CC
    del EE
    del multidarr

    # RRm = np.mean(RR,axis=2)
    # RRSTD = np.zeros_like(multidarr,dtype=np.float)
    # for mc in range(1, channels + 1):
    #     RRSTD[:,:,mc-1] = RR[:,:,mc-1] - RRm

    pvi = np.std(RR,axis=2)
    del RR
    # # print("max",frequency.max(),frequency.min())
    # return
    path = r'D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_Daily'

    path = os.path.join(path,"short_pvi_"+str(year)+".tif")
    pvi[mask>0] = -9999
    del mask
    write_Img(pvi,path,299,228)
    # return pvi
    del pvi


for i in range(2007,2008):
    PVI(r"D:\Cornell\EthiopianDrought\ChirpsDaily2", i)