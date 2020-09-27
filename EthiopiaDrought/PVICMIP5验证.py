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
                                   "cmip5_{}{}{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1
    print("bandNum",bandNum)
    mask = np.zeros((7, 9))
    multidarr = np.zeros((7, 9, bandNum))
    band_id = 0
    PList = []
    DVList = []

    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):
        # print(dt)
        chirps_file = os.path.join(chirpsdirectory,
                                   "cmip5_{}{}{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        # print(chirps_file)
        chirps = gdal.Open(chirps_file).ReadAsArray()
        PList.append(chirps[3,2])
        print(dt,chirps[3,2],"*")
        axisTime.append([dt])
        del chirps

    C = []
    temp = 0
    for i in range(len(PList)):
        temp += PList[i]
        C.append(temp)
    Pm = np.array(PList).mean()
    E = [(i+1)*Pm for i in range(len(PList))]

    R = [C[i]/E[i] for i in range(len(PList))]
    Rm = np.array(R).mean()

    Rc = [(R[i]-Rm)*(R[i]-Rm) for i in range(len(PList))]

    print(np.sqrt(np.array(Rc).mean()))




PVI(r"D:\Cornell\EthiopianDrought\CMIP5DailyClip", 2006)

#
# PList = []
# for i in range(60):
#     if i< 10:
#         print(i%6==0)
#         PList.append(10)
#     else:
#         PList.append(0)
# print(PList)
# C = []
# temp = 0
# for i in range(len(PList)):
#     temp += PList[i]
#     C.append(temp)
# Pm = np.array(PList).mean()
# E = [(i+1)*Pm for i in range(len(PList))]
#
# R = [C[i]/E[i] for i in range(len(PList))]
# Rm = np.array(R).mean()
#
# Rc = [(R[i]-Rm)*(R[i]-Rm) for i in range(len(PList))]
#
# print(np.sqrt(np.array(Rc).mean()))