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

def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float64):

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

def CMIP5(cmip,outputdir):
    start = datetime.strptime("-".join(["2006", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2020-12-31", "%Y-%m-%d").date()
    i = 0
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    geotrans = [0.9375,1.875,0,89.375,0,-1.25]
    for dt in (rrule.rrule(rrule.DAILY, interval=1, dtstart=start, until=stop)):
        data = cmip[i,1:144,:]
        data = np.flip(data,0)
        data2= data[:, 1:96]
        del data
        path = os.path.join(outputdir,"cmip5_{}{}{}.tif".format(str(dt.year),str(dt.month).zfill(2),str(dt.day).zfill(2)))
        write_Img(data2, path, proj, geotrans,95, 143,im_bands=1, dtype=gdal.GDT_Float64)
        del data2
        i += 1




#
path =r"C:\Users\zmhwh\Downloads\pr_day_ACCESS1-0_historicalExt_r2i1p1_20060101-20201231.nc"
#
fin = Dataset(path, "r")

pr = fin.variables["pr"][:]*3600*24
print(pr[0][2,2])
plt.imshow(pr[0])
plt.show()
# CMIP5(pr,r"D:\Cornell\EthiopianDrought\CMIP5Daily")


# data = pr.mean(axis=0)
#
# data = np.flip(data, 0)
# data = data[1:144, :]
#
#
# data2 = data[:, 1:96]
#
# proj = osr.SpatialReference()
# proj.ImportFromEPSG(4326)
# geotrans = [0.9375,1.875,0,89.375,0,-1.25]
# #
# path = os.path.join(r"C:\Users\zmhwh\Downloads","cmip5_mean2.tif")
# write_Img(data, path, proj, geotrans,192, 145,im_bands=1, dtype=gdal.GDT_Float32)
#
