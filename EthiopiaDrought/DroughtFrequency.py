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

def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

def AnomalyFrequency(chirpsdirectory,parameterspath,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    frequency = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))

    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        # print(chirps_file)
        paras = os.path.join(parameterspath, "chirps_month" + str(month).zfill(2) + ".tif")
        parameters = gdal.Open(paras)
        chirps = gdal.Open(chirps_file).ReadAsArray()

        averageMatrix = parameters.GetRasterBand(1).ReadAsArray()
        stdMatrix = parameters.GetRasterBand(2).ReadAsArray()
        maskarr = parameters.GetRasterBand(3).ReadAsArray()

        mask2 = (stdMatrix > 0) & (maskarr > -9999) & (chirps > -9999)
        mask2 = np.where(mask2)
        anomalyMap = np.zeros(shape=chirps.shape,dtype=np.float)
        anomalyMap[mask2] = (chirps[mask2] - averageMatrix[mask2]) / stdMatrix[mask2]
        frequency[anomalyMap < -1] +=1

        del chirps
        del parameters
        del averageMatrix
        del stdMatrix
        del maskarr
        del mask2
        del anomalyMap
    # print("max",frequency.max(),frequency.min())
    return frequency

