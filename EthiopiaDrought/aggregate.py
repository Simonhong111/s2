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
def aggregate500(input,output):

    stacks = glob.glob(os.path.join(input,'*.tif'))
    preference = gdal.Open(r'D:\Cornell\MCD12Q1v006\2003.tif')
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(preference.GetProjection()))
    geotrans = preference.GetGeoTransform()

    geo = [geotrans[0],geotrans[1]*2,geotrans[2],geotrans[3],geotrans[4],geotrans[5]*2]




    for stack in stacks:
        aggmatrix = np.full(shape=(2400, 2400), fill_value=255.0)
        dataset = gdal.Open(stack).ReadAsArray()
        for h in range(2400):
            for w in range(2400):
                if np.sum(dataset[(2*h):(2*h+2),(2*w):(2*w+2)] ==12) >2:
                    aggmatrix[h][w] = 12

        path = os.path.join(output, 'agg_' + os.path.basename(stack))
        write_Img(aggmatrix, path, proj, geo, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)



# aggregate500(r'D:\Cornell\MCD12Q1v006',r'D:\Cornell\MCD12Q1v006Agg')

def dayofyear(y,m,d):
    """
    :param y:
    :param m:
    :param d:
    :return: the nth day in the year
    """
    days_in_the_year = (date(y, m, d) - date(y, 1, 1)).days + 1
    return str(y)+str(days_in_the_year).zfill(3)
def aggr(directory,output):
    stacks = glob.glob(os.path.join(directory, '*.hdf'))

    for stack in stacks:
        aggmatrix = np.full(shape=(720, 1440), fill_value=255.0)
        ph21v07 = gdal.Open(stack).GetSubDatasets()[0][0]

        preference =gdal.Open(ph21v07)
        proj = osr.SpatialReference()
        proj.ImportFromWkt(str(preference.GetProjection()))
        geotrans = preference.GetGeoTransform()

        geo = [geotrans[0], geotrans[1] * 5, geotrans[2], geotrans[3], geotrans[4], geotrans[5] * 5]
        data = preference.ReadAsArray()
        for h in range(720):
            for w in range(1440):
                if np.sum(data[(5 * h):(5 * h + 5), (5 * w):(5 * w + 5)] == 12) > 12:
                    aggmatrix[h][w] = 12
        path = os.path.join(output, 'agg_' + os.path.basename(stack)+'.tif')
        write_Img(aggmatrix, path, proj, geo, 1440, 720, im_bands=1, dtype=gdal.GDT_Float32)

aggr(r"D:\Cornell\MCD12C1v006",r"D:\Cornell\MCD12C1v006")








    #
    # preference = gdal.Open(ph21v07)
    # proj = osr.SpatialReference()
    # proj.ImportFromWkt(str(preference.GetProjection()))
    # geotrans = preference.GetGeoTransform()
    # path = os.path.join(directory,str(dt.year)+'.tif')
    #
    # write_Img(MosaicImg, path, proj, geotrans, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)
    #
