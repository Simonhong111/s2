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
from epdatapro import write_Img,clipbyshp
def splitstack(input,output):

    stacks = glob.glob(os.path.join(input,'*'))

    preference = gdal.Open(r'D:\Cornell\NewSIF005Clip\SIF005_200208.nc.tif')
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(preference.GetProjection()))
    geotrans = preference.GetGeoTransform()
    W,H = preference.RasterXSize,preference.RasterYSize

    print(W,H)

    for stack in stacks:
        datasets = gdal.Open(stack)
        bandNum = datasets.RasterCount
        print('bandnum',bandNum)
        for i in range(bandNum):
            path = os.path.join(output, str(2007+i) + os.path.basename(stack))
            data = datasets.GetRasterBand(i+1).ReadAsArray()
            write_Img(data, path, proj, geotrans, W, H, im_bands=1, dtype=gdal.GDT_Float32)



splitstack(r'D:\Cornell\EthiopianDrought\AData\NewSIF0318Anom2007',r'D:\Cornell\EthiopianDrought\AData\NewSIF0318AnomSplit2007')



