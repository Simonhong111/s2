# This is used to generate seasonally Average
# of Precipitation (CHIRPS v2),EVI (MOD13C2),SIF (Jiaming)
# PVI (calculated from 10-day Precipitation)

import os,glob
import numpy as np
from dateutil import rrule
from datetime import *
from matplotlib import cm
from osgeo import gdal,osr,ogr
from matplotlib import pyplot as plt
from Experiment.Exp_Write_Image import write_Img as Write

def CHIRPSSeasonalComposite(ChirpsDir,OutDir,SeasonType="Short"):

    assert SeasonType in ["Short","Long"],"Season Type Error,select one from Short/Long"

    if SeasonType == "Short":
        MonthCollection = ["02","03","04","05"]
    if SeasonType == "Long":
        MonthCollection = ["06", "07", "08", "09"]

    MonthNumber = len(MonthCollection)

    RefPath = r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW,RefH = RefImg.RasterXSize,RefImg.RasterYSize

    # data cover 2003-2018
    for year in range(1981,2019):
        SeIMGs = np.zeros(shape=(RefH,RefW,MonthNumber),dtype=np.float)
        SeMask = np.zeros(shape=(RefH,RefW))

        for id,month in enumerate(MonthCollection):
            ch_path = os.path.join(ChirpsDir,"chirps-v2.0.{}.{}.tif".format(year,month))
            assert ch_path,"The {} does not exist".format(ch_path)
            print(ch_path)
            ch_data = gdal.Open(ch_path).ReadAsArray()
            SeIMGs[:,:,id] = ch_data
            SeMask[ch_data == -9999] = -9999

        SeAverage = SeIMGs.mean(axis=2)
        SeAverage[SeMask == -9999]

        out_path = os.path.join(OutDir,"chirps-v2.0.{}.{}.tif".format(year,SeasonType))
        Write(SeAverage, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)
        print("{} has been writen to disk".format(os.path.basename(out_path)))

CHIRPSSeasonalComposite(r"D:\Cornell\EthiopianDrought\Chirps2",
                        r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
                        SeasonType="Long")

def EVISeasonalComposite(EVIDir,OutDir,SeasonType="Short"):

    assert SeasonType in ["Short","Long"],"Season Type Error,select one from Short/Long"

    if SeasonType == "Short":
        MonthCollection = ["02","03","04","05"]
    if SeasonType == "Long":
        MonthCollection = ["06", "07", "08", "09"]

    MonthNumber = len(MonthCollection)

    RefPath = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW,RefH = RefImg.RasterXSize,RefImg.RasterYSize

    # data cover 2003-2018
    for year in range(2007,2019):
        SeIMGs = np.zeros(shape=(RefH,RefW,MonthNumber),dtype=np.float)
        SeMask = np.zeros(shape=(RefH,RefW))

        for id,month in enumerate(MonthCollection):
            Evi_path = os.path.join(EVIDir,"{}.{}.01.tif".format(year,month))
            assert Evi_path,"The {} does not exist".format(Evi_path)
            print(Evi_path)
            Evi_data = gdal.Open(Evi_path).ReadAsArray()
            SeIMGs[:,:,id] = Evi_data
            SeMask[Evi_data == -3000] = -9999
        SeAverage = SeIMGs.mean(axis=2)*(0.0001)
        SeAverage[SeMask == -9999] = -9999

        out_path = os.path.join(OutDir,"{}.{}.tif".format(year,SeasonType))
        Write(SeAverage, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)
        print("{} has been writen to disk".format(os.path.basename(out_path)))

# EVISeasonalComposite(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia",
#                         r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
#                         SeasonType="Short")

def SIFSeasonalComposite(SIFDir,OutDir,SeasonType="Short"):

    assert SeasonType in ["Short","Long"],"Season Type Error,select one from Short/Long"

    if SeasonType == "Short":
        MonthCollection = ["02","03","04","05"]
    if SeasonType == "Long":
        MonthCollection = ["06", "07", "08", "09"]

    MonthNumber = len(MonthCollection)

    # RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia\SIF005_200208.nc.tif"
    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia\SIF005_200208.nc.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW,RefH = RefImg.RasterXSize,RefImg.RasterYSize

    # data cover 2003-2018
    for year in range(2007,2019):
        SeIMGs = np.zeros(shape=(RefH,RefW,MonthNumber),dtype=np.float)
        SeMask = np.zeros(shape=(RefH,RefW))

        for id,month in enumerate(MonthCollection):
            SIF_path = os.path.join(SIFDir,"SIF005_{}{}.nc.tif".format(year,month))
            assert SIF_path,"The {} does not exist".format(SIF_path)
            print(SIF_path)
            SIF_data = gdal.Open(SIF_path).ReadAsArray()
            SeIMGs[:,:,id] = SIF_data
            SeMask[SIF_data < -0.05] = -9999
        SeAverage = SeIMGs.mean(axis=2)
        SeAverage[SeMask == -9999] = -9999

        out_path = os.path.join(OutDir,"SIF005_{}_{}.nc.tif".format(year,SeasonType))
        Write(SeAverage, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)
        print("{} has been writen to disk".format(os.path.basename(out_path)))

# SIFSeasonalComposite(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia",
#                         r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage",
#                         SeasonType="Long")

def PVISeasonalComposite(PVIDir,OutDir,SeasonType="short"):

    assert SeasonType in ["Short","Long"],"Season Type Error,select one from Short/Long"

    if SeasonType == "Short":
        MonthCollection = ["02","03","04","05"]
    if SeasonType == "Long":
        MonthCollection = ["06", "07", "08", "09"]

    MonthNumber = len(MonthCollection)

    # RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia\SIF005_200208.nc.tif"
    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_Month\pvi_200301.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW,RefH = RefImg.RasterXSize,RefImg.RasterYSize

    # data cover 2003-2018
    for year in range(2007,2019):
        SeIMGs = np.zeros(shape=(RefH,RefW,MonthNumber),dtype=np.float)
        SeMask = np.zeros(shape=(RefH,RefW))

        for id,month in enumerate(MonthCollection):
            PVI_path = os.path.join(PVIDir,"pvi_{}{}.tif".format(year,month))
            assert PVI_path,"The {} does not exist".format(PVI_path)
            print(PVI_path)
            PVI_data = gdal.Open(PVI_path).ReadAsArray()
            SeIMGs[:,:,id] = PVI_data
            SeMask[PVI_data == -9999] = -9999
        SeAverage = SeIMGs.mean(axis=2)
        SeAverage[SeMask == -9999] = -9999

        out_path = os.path.join(OutDir,"{}_pvi_{}.tif".format(SeasonType,year))
        Write(SeAverage, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)
        print("{} has been writen to disk".format(os.path.basename(out_path)))

# PVISeasonalComposite(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_Month",
#                         r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_Season",
#                         SeasonType="Long")


