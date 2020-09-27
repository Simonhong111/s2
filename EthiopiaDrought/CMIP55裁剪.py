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

def chirpsclip(chirpsdirectory,clippeddirectory,dstNodata=-9999):


    ref = gdal.Open(r"D:\Cornell\EthiopianDrought\CMIP5Daily\cmip5_20060101.tif")
    geo = ref.GetGeoTransform()
    proj = ref.GetProjection()
    ulx = geo[0] + 16*geo[1]
    uly = geo[3] + 58*geo[5]
    # 30 70
    # 14
    # 12
    geotrans =[ulx,1.875,0,uly,0,-1.25]
    chirps_decompress_files = glob.glob(os.path.join(chirpsdirectory,"*.tif"))
    # print(chirps_decompress_files)
    for tiffile in chirps_decompress_files:
        raster = gdal.Open(tiffile).ReadAsArray(16,58,10,12)
        path = os.path.join(clippeddirectory,os.path.basename(tiffile))
        print(path)
        write_Img(raster, path, proj, geotrans, 10, 12, im_bands=1, dtype=gdal.GDT_Float32)


chirpsclip(r"D:\Cornell\EthiopianDrought\CMIP5Daily",
              "D:\Cornell\EthiopianDrought\CMIP5DailyClip")




