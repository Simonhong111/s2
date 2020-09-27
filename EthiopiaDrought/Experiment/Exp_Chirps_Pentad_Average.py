from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
from calendar import monthrange

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


def PentandAverage(pentandDir,outputdirectory):
    refpath = r"D:\Cornell\EthiopianDrought\CHIRPS5Clip\chirps-v2.0.2003.01.1.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    With,Height = reference.RasterXSize,reference.RasterYSize
    chirp_files = glob.glob(r"D:\Cornell\EthiopianDrought\CHIRPS5Clip\chirps-v2.0.*.tif")

    for file in chirp_files:
        Month = os.path.basename(file)[-8:-6]
        Year = os.path.basename(file)[-13:-9]
        timeT = datetime.strptime("-".join([Year, Month, "01"]), "%Y-%m-%d")
        Days = monthrange(timeT.year,timeT.month)
        RestDay = Days[1] - 25
        chirp_data = gdal.Open(file).ReadAsArray()
        mask = np.where(chirp_data == -9999)
        if os.path.basename(file)[-5] =="6":
            chirp_data = chirp_data/RestDay

        if os.path.basename(file)[-5] in ["1","2","3","4","5"]:
            chirp_data = chirp_data / 5.0

        chirp_data[mask] = -9999


        path = os.path.join(outputdirectory, os.path.basename(file))
        write_Img(chirp_data, path, proj, geotrans, With, Height, im_bands=1, dtype=gdal.GDT_Float32)


PentandAverage(r"D:\Cornell\EthiopianDrought\CHIRPS5Clip",
               r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_Pentad_Average")



