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

def chirpsclip(chirpsdirectory,epregionshppath,clippeddirectory,dstNodata=-9999):

    chirps_decompress_files = glob.glob(os.path.join(chirpsdirectory,"*.tif"))
    print(chirps_decompress_files)
    for tiffile in chirps_decompress_files:
        input_raster = tiffile
        output_raster = os.path.join(clippeddirectory,os.path.basename(input_raster))
        clipbyshp(input_raster,output_raster,epregionshppath,dstNodata=dstNodata)
        print("{} has been processed".format(input_raster))


chirpsclip(r"D:\Cornell\EthiopianDrought\CHIRPS5",
              r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
              "D:\Cornell\EthiopianDrought\CHIRPS5Clip")



# start = datetime.strptime("-".join(["2003", "01", "01"]), "%Y-%m-%d").date()
# stop = datetime.strptime("-".join(["2018", "12", "31"]), "%Y-%m-%d").date()

# baseyear = 2003
# for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
#     for day in range(1,4):
#         path = os.path.join("D:\Cornell\EthiopianDrought\CHIRPS10",
#                             "chirps-v2.0.{}.{}.{}.tif".format(str(dt.year), str(dt.month).zfill(2),
#                                                               str(day)))
#         if os.path.exists(path):
#             print(path)


