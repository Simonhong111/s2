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
def clipbyshp(input_raster,output_raster,input_shape, dstNodata=-9999):
    """
    :param input_raster: the raster data being processed later
    :param output_raster: the clipped datas' savepaths
    :param input_shape: the shape defining the extent
    :return: none
    """
    ds = gdal.Warp(output_raster,
                   input_raster,
                   format='GTiff',
                   cutlineDSName=input_shape,  # or any other file format
                   # cutlineDSName=None,
                   # cutlineWhere="FIELD = 'whatever'",
                   # optionally you can filter your cutline (shapefile) based on attribute values
                   cropToCutline=True,
                   dstNodata=dstNodata)  # select the no data value you like
    ds = None

def write_Img(data, path, im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float64):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    rainfall = gdal.Open(r"D:\Cornell\EthiopianDrought\CMIP5DailyClip\cmip5_20060101.tif")
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
                                   "cmip5_{}{}{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1
    print("bandNum",bandNum)
    mask = np.zeros((12, 10))
    multidarr = np.zeros((12, 10, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):
        print(dt)
        chirps_file = os.path.join(chirpsdirectory,
                                   "cmip5_{}{}{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        print(chirps_file)
        chirps = gdal.Open(chirps_file).ReadAsArray()

        multidarr[:,:,band_id] = chirps
        mask[chirps == -9999] += 1
        mask[chirps >10000] +=1
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

    # RRm = np.mean(RR,axis=2)
    # RRSTD = np.zeros_like(multidarr,dtype=np.float)
    # for mc in range(1, channels + 1):
    #     RRSTD[:,:,mc-1] = RR[:,:,mc-1] - RRm

    pvi = np.std(RR,axis=2)

    # # print("max",frequency.max(),frequency.min())
    # return
    path = r'D:\Cornell\EthiopianDrought\AData\CMIP5PVI'

    path = os.path.join(path,"long_pvi_"+str(year)+".tif")
    pvi[ppm==0] = 0
    pvi[mask>0] = -9999
    write_Img(pvi,path,10,12)
    # return pvi

for i in range(2006,2019):
    PVI(r"D:\Cornell\EthiopianDrought\CMIP5DailyClip", i)