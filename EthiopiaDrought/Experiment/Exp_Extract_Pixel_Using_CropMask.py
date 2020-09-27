from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt


def ExtractPixels(ChirpsDir,EviDir,SIFDir,PVIDir,CropPath,StartYear,EndYear):
    """
    :param ChirpsDir: precipitation directory
    :param EviDir:
    :param SIFDir:
    :param PVIDir:
    :param CropPath: crop mask
    :return:
    """
    crop = gdal.Open(CropPath)
    Width,Height = crop.RasterXSize,crop.RasterYSize
    crop_raster = crop.ReadAsArray()
    mask = np.where(crop_raster == 255)  # get crop mask locations

    vRow = mask[0]
    vCol = mask[1]
    ind = vRow*Width + vCol  # the location of crop pixel in the crop mask
    RF_SHORT = []
    RF_LONG = []
    EVI_SHORT = []
    EVI_LONG = []
    SIF_SHORT = []
    SIF_LONG = []
    PVI_SHORT = []
    PVI_LONG = []
    COLSet = []
    ROWSet = []

    for year in range(StartYear,EndYear+1):
        rf_file_short = os.path.join(ChirpsDir,"chirps-v2.0.{}.Short.tif".format(year))
        rf_file_long = os.path.join(ChirpsDir, "chirps-v2.0.{}.Long.tif".format(year))
        evi_file_short = os.path.join(EviDir,"{}.Short.tif".format(year))
        evi_file_long = os.path.join(EviDir, "{}.Long.tif".format(year))
        sif_file_short = os.path.join(SIFDir,"SIF005_{}_Short.nc.tif".format(year))
        sif_file_long = os.path.join(SIFDir, "SIF005_{}_Long.nc.tif".format(year))
        pvi_file_short = os.path.join(PVIDir, "short_pvi_{}.tif".format(year))
        pvi_file_long = os.path.join(PVIDir, "long_pvi_{}.tif".format(year))

        assert os.path.exists(rf_file_short),"the file does not exists {}".format(rf_file_short)
        assert os.path.exists(rf_file_long), "the file does not exists {}".format(rf_file_long)
        assert os.path.exists(evi_file_short), "the file does not exists {}".format(evi_file_short)
        assert os.path.exists(evi_file_long), "the file does not exists {}".format(evi_file_long)
        assert os.path.exists(sif_file_short), "the file does not exists {}".format(sif_file_short)
        assert os.path.exists(sif_file_long), "the file does not exists {}".format(sif_file_long)
        assert os.path.exists(pvi_file_short), "the file does not exists {}".format(pvi_file_short)
        assert os.path.exists(pvi_file_long), "the file does not exists {}".format(pvi_file_long)


        rf_data_short = gdal.Open(rf_file_short).ReadAsArray()
        rf_data_long = gdal.Open(rf_file_long).ReadAsArray()
        evi_data_short = gdal.Open(evi_file_short).ReadAsArray()
        evi_data_long = gdal.Open(evi_file_long).ReadAsArray()
        sif_data_short = gdal.Open(sif_file_short).ReadAsArray()
        sif_data_long = gdal.Open(sif_file_long).ReadAsArray()
        pvi_data_short = gdal.Open(pvi_file_short).ReadAsArray()
        pvi_data_long = gdal.Open(pvi_file_long).ReadAsArray()


        rf_crop_short = np.take(rf_data_short,ind)
        rf_crop_long = np.take(rf_data_long, ind)
        evi_crop_short = np.take(evi_data_short, ind)
        evi_crop_long = np.take(evi_data_long, ind)
        sif_crop_short = np.take(sif_data_short, ind)
        sif_crop_long = np.take(sif_data_long, ind)
        pvi_crop_short = np.take(pvi_data_short, ind)
        pvi_crop_long = np.take(pvi_data_long, ind)

        RF_SHORT.append(rf_crop_short.tolist())
        RF_LONG.append(rf_crop_long.tolist())
        EVI_SHORT.append(evi_crop_short.tolist())
        EVI_LONG.append(evi_crop_long.tolist())
        SIF_SHORT.append(sif_crop_short.tolist())
        SIF_LONG.append(sif_crop_long.tolist())
        PVI_SHORT.append(pvi_crop_short.tolist())
        PVI_LONG.append(pvi_crop_long.tolist())

        del rf_data_short
        del rf_data_long
        del evi_data_short
        del evi_data_long
        del sif_data_short
        del sif_data_long
        del pvi_data_short
        del pvi_data_long

        del rf_crop_short
        del rf_crop_long
        del evi_crop_short
        del evi_crop_long
        del sif_crop_short
        del sif_crop_long
        del pvi_crop_short
        del pvi_crop_long

    COLSet.append(vCol)
    ROWSet.append(vRow)

    return RF_SHORT,RF_LONG,EVI_SHORT,EVI_LONG,SIF_SHORT,SIF_LONG,PVI_SHORT,PVI_LONG,COLSet,ROWSet


ChirpsDir = r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage"
EviDir = r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage"
SIFDir = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage"
PVIDir = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D"
CropPath  = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropMask\agg_clip50.tif"
StartYear = 2007
EndYear = 2018
RF_SHORT,RF_LONG,EVI_SHORT,EVI_LONG,SIF_SHORT,SIF_LONG,PVI_SHORT,PVI_LONG,COLSet,ROWSet = ExtractPixels(ChirpsDir,EviDir,SIFDir,PVIDir,CropPath,StartYear,EndYear)
YearNum = EndYear - StartYear + 1
RFD_S = np.array(RF_SHORT).T.reshape(-1,YearNum)
RFD_L = np.array(RF_LONG).T.reshape(-1,YearNum)
EVID_S = np.array(EVI_SHORT).T.reshape(-1,YearNum)
EVID_L = np.array(EVI_LONG).T.reshape(-1,YearNum)
SIFD_S = np.array(SIF_SHORT).T.reshape(-1,YearNum)
SIFD_L = np.array(SIF_LONG).T.reshape(-1,YearNum)
PVID_S = np.array(PVI_SHORT).T.reshape(-1,YearNum)
PVID_L = np.array(PVI_LONG).T.reshape(-1,YearNum)
ROWSetList = np.array(ROWSet).T.reshape(-1,1)
COLSetList = np.array(COLSet).T.reshape(-1,1)


mask = np.zeros(SIFD_S.shape[0],dtype=np.int)
for id in range(SIFD_S.shape[1]):
    mask[np.where(SIFD_S[:,id]<-0.05)] = -9999
    mask[np.where(SIFD_L[:,id]<-0.05)] = -9999
valid = np.where(mask !=-9999) # delete the value  less than -0.05

print(len(valid[0]))



DataDict = {}

for index,year in enumerate(range(StartYear,EndYear+1)):
    DataDict["RF" + str(year) + "Short"] = RFD_S[:, index][valid]
    DataDict["RF" + str(year) + "Long"] = RFD_L[:, index][valid]
    DataDict["EVI"+str(year)+"Short"] = EVID_S[:,index][valid]
    DataDict["EVI" + str(year) + "Long"] = EVID_L[:, index][valid]
    DataDict["SIF" + str(year) + "Short"] = SIFD_S[:, index][valid]
    DataDict["SIF" + str(year) + "Long"] = SIFD_L[:, index][valid]
    DataDict["PVI" + str(year) + "Short"] = PVID_S[:, index][valid]
    DataDict["PVI" + str(year) + "Long"] = PVID_L[:, index][valid]

DataDict["ROW"] = ROWSetList[valid].flatten()
DataDict["COL"] = COLSetList[valid].flatten()


df = pd.DataFrame(DataDict)


outpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"
df.to_csv(outpath, index=False)


