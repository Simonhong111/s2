from osgeo import gdal,osr,ogr
import os
import numpy as np
from dateutil import rrule
from datetime import *
from scipy import signal
from Experiment.Exp_Write_Image import write_Img
import ctypes
ctypes.WinDLL("kernel32.dll")
ctypes.WinDLL("msvcrt.dll")
ctypes.WinDLL("user32.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libwinpthread-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgcc_s_seh-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgomp-1.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\cblas\.libs\libgslcblas-0.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\.libs\libgsl-25.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\libeemd\libeemd.so")
from pyeemd import ceemdan,eemd
from pyeemd.utils import plot_imfs
import matplotlib.pyplot as plt

def fit(y):
    # imfs = ceemdan(y,num_imfs = 4,S_number=4,num_siftings=50)
    # trend = imfs[-1]
    # dtrend = y-trend

    return signal.detrend(y)

def detrendSIF(SIFDir,OutDir,SeasonType):
    st_year = 2007
    refpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage\SIF005_2007_Long.nc.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    Width, Height = reference.RasterXSize, reference.RasterYSize

    bandNum = 0

    for year in range(st_year,2019):

        sif_file = os.path.join(SIFDir,
                                   "SIF005_{}_{}.nc.tif".format(str(year), SeasonType))
        if os.path.exists(sif_file):
            bandNum += 1
        assert os.path.exists(sif_file),"{} does not exist".format(sif_file)

    maskarr = np.zeros((Height, Width))
    multidarr = np.zeros((Height, Width, bandNum),dtype=np.float)
    band_id = 0

    for year in range(st_year,2019):

        sif_file = os.path.join(SIFDir,
                                   "SIF005_{}_{}.nc.tif".format(year, SeasonType))

        if os.path.exists(sif_file):
            sif_data = gdal.Open(sif_file).ReadAsArray()
            maskarr[sif_data < -0.05] = -9999
            multidarr[:, :, band_id] = sif_data
            band_id += 1

    multiarrflatten = multidarr.reshape(Width*Height,bandNum)
    del multidarr
    print("*")
    results = map(fit,multiarrflatten)
    del multiarrflatten
    print("**")
    deResults = np.array(list(results)).reshape(Height,Width,bandNum)
    del results

    for i in range(bandNum):
        deResults[:,:,i][maskarr ==-9999] = -9999
        outputpath = os.path.join(OutDir, "SIF005_{}_{}.nc.tif".format(i +st_year,SeasonType))
        write_Img(deResults[:,:,i], outputpath, proj, geotrans, Width, Height, im_bands=1, dtype=gdal.GDT_Float32)


# res = detrendSIF(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage",
#                      r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage","Short")
from multiprocessing import Pool
if __name__ == '__main__':
    st_year = 2007
    refpath =r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage\SIF005_2007_Long.nc.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    Width, Height = reference.RasterXSize, reference.RasterYSize
    SIFDir =r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage"
    OutDir = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage"
    SeasonType = "Short"


    bandNum = 0

    for year in range(st_year, 2019):

        sif_file = os.path.join(SIFDir,
                                "SIF005_{}_{}.nc.tif".format(str(year), SeasonType))

        if os.path.exists(sif_file):
            bandNum += 1
        assert os.path.exists(sif_file), "{} does not exist".format(sif_file)

    maskarr = np.zeros((Height, Width))
    multidarr = np.zeros((Height, Width, bandNum), dtype=np.float)
    band_id = 0

    for year in range(st_year, 2019):

        sif_file = os.path.join(SIFDir,
                                "SIF005_{}_{}.nc.tif".format(year, SeasonType))

        if os.path.exists(sif_file):
            sif_data = gdal.Open(sif_file).ReadAsArray()
            maskarr[sif_data < -0.05] = -9999
            multidarr[:, :, band_id] = sif_data
            band_id += 1

    multiarrflatten = multidarr.reshape(Width * Height, bandNum)
    del multidarr
    print("*")
    with Pool(5) as p:
        results = p.map(fit, multiarrflatten)
    del multiarrflatten
    print("**")
    deResults = np.array(list(results)).reshape(Height, Width, bandNum)
    del results

    for i in range(bandNum):
        deResults[:, :, i][maskarr == -9999] = -9999
        outputpath = os.path.join(OutDir, "SIF005_{}_{}.nc.tif".format(i + st_year, SeasonType))
        write_Img(deResults[:, :, i], outputpath, proj, geotrans, Width, Height, im_bands=1, dtype=gdal.GDT_Float32)


