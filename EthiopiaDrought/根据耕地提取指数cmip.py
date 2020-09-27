from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt
# ET\

def ExtractRainFall():

    start = datetime.strptime("-".join(["2006", "01", "01"]), "%Y-%m-%d").date()
    stop =  datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    crop_path = r"D:\Cornell\EthiopianDrought\CropType2015\cmip_agg_clip50.tif"
    crop = gdal.Open(crop_path)
    Width,Height = crop.RasterXSize,crop.RasterYSize
    crop_raster = crop.ReadAsArray()
    mask = np.where(crop_raster == 255)

    vRow = mask[0]
    vCol = mask[1]
    ind = vRow*Width + vCol

    CMIPRF_Short = []
    CMIPRF_Long = []
    CMIPPVI_Short = []
    CMIPPVI_Long = []

    CHIRPRF_Short = []
    CHIRPRF_Long = []
    CHIRPPVI_Short = []
    CHIRPPVI_Long = []

    COLSet = []
    ROWSet = []

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        short_cmip_rf = os.path.join(r"D:\Cornell\EthiopianDrought\CMIPMonth",
                                  "cmip5_short_{}.tif".format(dt.year))
        long_cmip_rf = os.path.join(r"D:\Cornell\EthiopianDrought\CMIPMonth",
                                     "cmip5_long_{}.tif".format(dt.year))
        short_cmip_pvi = os.path.join(r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI",
                                     "short_pvi_{}.tif".format(dt.year))
        long_cmip_pvi = os.path.join(r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI",
                                     "long_pvi_{}.tif".format(dt.year))

        short_chirp_rf = os.path.join(r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth",
                                     "chirps-v2.0.short_{}.tif".format(dt.year))
        long_chirp_rf = os.path.join(r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth",
                                    "chirps-v2.0.long_{}.tif".format(dt.year))
        short_chirp_pvi = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVIDaily",
                                      "short_pvi_{}.tif".format(dt.year))
        long_chirp_pvi = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVIDaily",
                                     "long_pvi_{}.tif".format(dt.year))


        assert os.path.exists(short_cmip_rf), "the file does not exists {}".format(short_cmip_rf)
        assert os.path.exists(long_cmip_rf), "the file does not exists {}".format(long_cmip_rf)
        assert os.path.exists(short_cmip_pvi), "the file does not exists {}".format(short_cmip_pvi)
        assert os.path.exists(long_cmip_pvi), "the file does not exists {}".format(long_cmip_pvi)

        assert os.path.exists(short_chirp_rf), "the file does not exists {}".format(short_chirp_rf)
        assert os.path.exists(long_chirp_rf), "the file does not exists {}".format(long_chirp_rf)
        assert os.path.exists(short_chirp_pvi), "the file does not exists {}".format(short_chirp_pvi)
        assert os.path.exists(long_chirp_pvi), "the file does not exists {}".format(long_chirp_pvi)

        short_cmip_rf_data = gdal.Open(short_cmip_rf).ReadAsArray()
        long_cmip_rf_data = gdal.Open(long_cmip_rf).ReadAsArray()
        short_cmip_pvi_data = gdal.Open(short_cmip_pvi).ReadAsArray()
        long_cmip_pvi_data = gdal.Open(long_cmip_pvi).ReadAsArray()

        short_chirp_rf_data = gdal.Open(short_chirp_rf).ReadAsArray()
        long_chirp_rf_data = gdal.Open(long_chirp_rf).ReadAsArray()
        short_chirp_pvi_data = gdal.Open(short_chirp_pvi).ReadAsArray()
        long_chirp_pvi_data = gdal.Open(long_chirp_pvi).ReadAsArray()

        short_cmip_rf_V = np.take(short_cmip_rf_data, ind)
        long_cmip_rf_V = np.take(long_cmip_rf_data, ind)
        short_cmip_pvi_V = np.take(short_cmip_pvi_data, ind)
        long_cmip_pvi_V = np.take(long_cmip_pvi_data, ind)

        short_chirp_rf_V = np.take(short_chirp_rf_data, ind)
        long_chirp_rf_V = np.take(long_chirp_rf_data, ind)
        short_chirp_pvi_V = np.take(short_chirp_pvi_data, ind)
        long_chirp_pvi_V = np.take(long_chirp_pvi_data, ind)

        CMIPRF_Short.append(short_cmip_rf_V )
        CMIPRF_Long.append(long_cmip_rf_V)
        CMIPPVI_Short.append(short_cmip_pvi_V)
        CMIPPVI_Long.append(long_cmip_pvi_V)

        CHIRPRF_Short.append(short_chirp_rf_V )
        CHIRPRF_Long.append(long_chirp_rf_V)
        CHIRPPVI_Short.append(short_chirp_pvi_V)
        CHIRPPVI_Long.append(long_chirp_pvi_V)


    COLSet.append(vCol)
    ROWSet.append(vRow)

    return CMIPRF_Short,CMIPRF_Long,CMIPPVI_Short,CMIPPVI_Long,CHIRPRF_Short,CHIRPRF_Long,CHIRPPVI_Short,CHIRPPVI_Long,COLSet,ROWSet

CMIPRF_Short,CMIPRF_Long,CMIPPVI_Short,CMIPPVI_Long,CHIRPRF_Short,CHIRPRF_Long,CHIRPPVI_Short,CHIRPPVI_Long,COLSet,ROWSet = ExtractRainFall()


CMIPRF_Short = np.array(CMIPRF_Short).T.reshape(-1,13)
CMIPRF_Long = np.array(CMIPRF_Long).T.reshape(-1,13)
CMIPPVI_Short = np.array(CMIPPVI_Short).T.reshape(-1,13)
CMIPPVI_Long = np.array(CMIPPVI_Long ).T.reshape(-1,13)

CHIRPRF_Short = np.array(CHIRPRF_Short).T.reshape(-1,13)
CHIRPRF_Long = np.array(CHIRPRF_Long).T.reshape(-1,13)
CHIRPPVI_Short = np.array(CHIRPPVI_Short).T.reshape(-1,13)
CHIRPPVI_Long = np.array(CHIRPPVI_Long ).T.reshape(-1,13)


ROWSetList = np.array(ROWSet).T.reshape(-1,1)
COLSetList = np.array(COLSet).T.reshape(-1,1)



Month =[1,2,3,4,5,6,7,8,9,10,11,12]
Year =range(2006,2019)

DataDict = {}



for index,year in enumerate(Year):
    DataDict["ShortCMRF" +str(year)] = CMIPRF_Short[:, index]
    DataDict["LongCMRF" + str(year)] = CMIPRF_Long[:, index]
    DataDict["ShortCMPVI" + str(year)] = CMIPPVI_Short[:, index]
    DataDict["LongCMPVI" + str(year)] = CMIPPVI_Long[:, index]

    DataDict["ShortCHRF" + str(year)] = CHIRPRF_Short[:, index]
    DataDict["LongCHRF" + str(year)] = CHIRPRF_Long[:, index]
    DataDict["ShortCHPVI" + str(year)] = CHIRPPVI_Short[:, index]
    DataDict["LongCHPVI" + str(year)] = CHIRPPVI_Long[:, index]

DataDict["ROW"] = ROWSetList.flatten()
DataDict["COL"] = COLSetList.flatten()



df = pd.DataFrame(DataDict)
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\CMIP\CMIP_PolyGonAgg_Mask50.csv"
df.to_csv(outpath, index=False)

