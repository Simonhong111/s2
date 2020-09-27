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
    path = r"D:\Cornell\EthiopianDrought\Chirps2"
    start = datetime.strptime("-".join(["2003", "01", "01"]), "%Y-%m-%d").date()
    stop =  datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    crop_path = r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip70.tif"
    crop = gdal.Open(crop_path)
    Width,Height = crop.RasterXSize,crop.RasterYSize
    crop_raster = crop.ReadAsArray()
    mask = np.where(crop_raster == 255)
    # print(mask)
    vRow = mask[0]
    vCol = mask[1]
    ind = vRow*Width + vCol
    EVI = []
    RF = []
    NSIF = []
    GSIF = []
    COLSet = []
    ROWSet = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        rf_file = os.path.join(path,"chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        gsif_file = os.path.join(r"D:\Cornell\GOSIFV002Clip","GOSIF_{}.M{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        nsif_file = os.path.join(r"D:\Cornell\NewSIF005Clip","SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))
        evi_file = os.path.join(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia","{}.{}.01.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(rf_file),"the file does not exists {}".format(rf_file)
        assert os.path.exists(gsif_file), "the file does not exists {}".format(gsif_file)
        assert os.path.exists(nsif_file), "the file does not exists {}".format(nsif_file)
        assert os.path.exists(evi_file), "the file does not exists {}".format(evi_file)
        # print(chirps_file)
        RFData = gdal.Open(rf_file).ReadAsArray()
        GSIFData = gdal.Open(gsif_file).ReadAsArray()
        NSIFData = gdal.Open(nsif_file).ReadAsArray()
        EVIData = gdal.Open(evi_file).ReadAsArray()

        RFV = np.take(RFData,ind)
        GSIFV = np.take(GSIFData, ind)
        NSIFV = np.take(NSIFData, ind)
        EVIV = np.take(EVIData, ind)

        EVI.append(EVIV.tolist())
        RF.append(RFV.tolist())
        NSIF.append(NSIFV.tolist())
        GSIF.append(GSIFV.tolist())
        del RFData
        del GSIFData
        del NSIFData
        del EVIData
        del RFV
        del GSIFV
        del NSIFV
        del EVIV
    COLSet.append(vCol)
    ROWSet.append(vRow)

    SPVI = []
    LPVI = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        short_file = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVI5day",
                                  "short_pvi_{}.tif".format(dt.year))
        long_file = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVI5day",
                                 "long_pvi_{}.tif".format(dt.year))

        assert os.path.exists(short_file), "the file does not exists {}".format(short_file)
        assert os.path.exists(long_file), "the file does not exists {}".format(long_file)


        ShortData = gdal.Open(short_file).ReadAsArray()
        LongData = gdal.Open(long_file).ReadAsArray()


        ShortV = np.take(ShortData, ind)
        LongV = np.take(LongData, ind)


        SPVI.append(ShortV.tolist())
        LPVI.append(LongV.tolist())

        del ShortData
        del LongData


    return EVI,RF,NSIF,GSIF,SPVI,LPVI,COLSet,ROWSet

EVI, RF, NSIF, GSIF, SPVI, LPVI, COLSet, ROWSet = ExtractRainFall()


EVIDList = np.array(EVI).T.reshape(-1,192)
RFDList = np.array(RF).T.reshape(-1,192)
NSIFDList = np.array(NSIF).T.reshape(-1,192)
GSIFDList = np.array(GSIF).T.reshape(-1,192)
SPVIDList = np.array(SPVI).T.reshape(-1,16)
LPVIDList = np.array(LPVI).T.reshape(-1,16)
ROWSetList = np.array(ROWSet).T.reshape(-1,1)
COLSetList = np.array(COLSet).T.reshape(-1,1)

evimask = np.where(EVIDList == -3000)
EVIDList = EVIDList*0.0001
EVIDList[evimask] = -9999
gsifmask = np.where(GSIFDList >=32766)
GSIFDList = GSIFDList*0.0001
GSIFDList[gsifmask] = -9999

print(NSIFDList.shape)
mask = np.zeros(NSIFDList.shape[0],dtype=np.int)
for id in range(NSIFDList.shape[1]):
    mask[np.where(NSIFDList[:,id]<-0.05)] = 1
valid = np.where(mask !=1)




Month =[1,2,3,4,5,6,7,8,9,10,11,12]
Year =range(2003,2019)
YM = []
for year in Year:
    for m in Month:
        YM.append(str(year)+str(m).zfill(2))
DataDict = {}

for index,ym in enumerate(YM):
    DataDict["EVI"+ym] = EVIDList[:,index][valid]
    DataDict["RF" + ym] = RFDList[:, index][valid]
    DataDict["NSIF" + ym] = NSIFDList[:, index][valid]
    DataDict["GSIF" + ym] = GSIFDList[:, index][valid]


for index,year in enumerate(Year):
    DataDict["ShortPVI" +str(year)] = SPVIDList[:, index][valid]
    DataDict["LongPVI" + str(year)] = LPVIDList[:, index][valid]
DataDict["ROW"] = ROWSetList[valid].flatten()
DataDict["COL"] = COLSetList[valid].flatten()



df = pd.DataFrame(DataDict)
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\5day\PolyGonAgg_Mask70.csv"
df.to_csv(outpath, index=False)

