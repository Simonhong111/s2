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

    rainfall = gdal.Open(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_Pentad_Aggragation\agg_chirps-v2.0.2003.01.1.tif")
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

def PVI(cmip5directory,year):
    start = datetime.strptime("-".join([str(year), str(6).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("-".join([str(year), str(9).zfill(2), "30"]), "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):

        for day in range(1,7):
            chirps_file = os.path.join(cmip5directory,
                                       "agg_chirps-v2.0.{}.{}.{}.tif".format(str(dt.year), str(dt.month).zfill(2),
                                                                         str(day)))
            if os.path.exists(chirps_file):
                bandNum += 1
    print("bandNum",bandNum)
    mask = np.zeros((20, 13))
    multidarr = np.zeros((20, 13, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):

        for day in range(1,7):
            chirps_file = os.path.join(cmip5directory,
                                       "agg_chirps-v2.0.{}.{}.{}.tif".format(str(dt.year), str(dt.month).zfill(2),
                                                                         str(day)))
            print(chirps_file)
            chirps = gdal.Open(chirps_file).ReadAsArray()

            multidarr[:, :, band_id] = chirps
            mask[chirps == -9999] += 1
            band_id += 1
            axisTime.append([dt.year, dt.month])

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

    # RRm = np.mean(RR,axis=2)
    # RRSTD = np.zeros_like(multidarr,dtype=np.float)
    # for mc in range(1, channels + 1):
    #     RRSTD[:,:,mc-1] = RR[:,:,mc-1] - RRm

    pvi = np.std(RR,axis=2)

    # # print("max",frequency.max(),frequency.min())
    # return
    path = r'D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CHIRPS'

    path = os.path.join(path,"long_pvi_"+str(year)+".tif")
    pvi[mask>0] = -9999
    write_Img(pvi,path,13,20)
    # return pvi

for i in range(2007,2019):
    PVI(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_Pentad_Aggragation", i)