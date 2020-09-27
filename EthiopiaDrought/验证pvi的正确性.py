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
    start = datetime.strptime("-".join([str(year), str(6).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("-".join([str(year), str(9).zfill(2), "30"]), "%Y-%m-%d").date()
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
    PList = []
    DVList = []

    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):
        # print(dt)
        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0.{}.{}.{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        # print(chirps_file)
        chirps = gdal.Open(chirps_file).ReadAsArray()
        PList.append(chirps[202,179])
        print(chirps[202,179],"*")
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




# PVI(r"D:\Cornell\EthiopianDrought\ChirpsDaily2", 2003)


PList = []
for i in range(60):
    if i< 56:
        print(i%6==0)
        PList.append(0)
    else:
        PList.append(15)
print(PList)
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