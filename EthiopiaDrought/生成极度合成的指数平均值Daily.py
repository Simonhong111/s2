from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
def monthAnalysis(path,monthtype,var1,var2,year):
    if monthtype == "Short":
        months = [2,3,4,5]
    if monthtype == "Long":
        months =[6,7,8,9]
    data = pd.read_csv(path)
    key1 = [var1+str(year)+str(months[0]).zfill(2),var1+str(year)+str(months[1]).zfill(2),
            var1+str(year)+str(months[2]).zfill(2),var1+str(year)+str(months[3]).zfill(2)]
    key2 = [var2 + str(year) + str(months[0]).zfill(2), var2 + str(year) + str(months[1]).zfill(2),
            var2 + str(year) + str(months[2]).zfill(2), var2 + str(year) + str(months[3]).zfill(2)]
    # print(data[key1].to_numpy())
    # print(data[key1].to_numpy().mean(axis=1))
    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    return Var1List,Var2List,PVIList



path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask70.csv"
Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "EVI", "RF", 2010)
row_num = Var1List.shape[0]
col_num = 16

evi_short =np.zeros(shape=(row_num,col_num),dtype=np.float)
rf_short = np.zeros(shape=(row_num,col_num),dtype=np.float)
nsif_short = np.zeros(shape=(row_num,col_num),dtype=np.float)
gsif_short = np.zeros(shape=(row_num,col_num),dtype=np.float)
spvi = np.zeros(shape=(row_num,col_num),dtype=np.float)


evi_long = np.zeros(shape=(row_num,col_num),dtype=np.float)
rf_long = np.zeros(shape=(row_num,col_num),dtype=np.float)
nsif_long = np.zeros(shape=(row_num,col_num),dtype=np.float)
gsif_long = np.zeros(shape=(row_num,col_num),dtype=np.float)
lpvi = np.zeros(shape=(row_num,col_num),dtype=np.float)
for id,Year in enumerate(range(2003,2019)):
    Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "EVI", "RF", Year)
    evi_short[:,id] = Var1List
    rf_short[:,id] = Var2List
    spvi[:,id] = PVIList

    Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "EVI", "RF", Year)
    evi_long[:, id] = Var1List
    rf_long[:, id] = Var2List
    lpvi[:, id] = PVIList

    Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "GSIF", "RF", Year)
    gsif_short[:, id] = Var1List

    Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "GSIF", "RF", Year)
    gsif_long[:, id] = Var1List

    Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "NSIF", "RF", Year)
    nsif_short[:, id] = Var1List

    Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "NSIF", "RF", Year)
    nsif_long[:, id] = Var1List


#
rf_short2 = (rf_short - rf_short.mean(axis=1)[:,np.newaxis])/rf_short.std(axis=1)[:,np.newaxis]
rf_long2 = (rf_long - rf_long.mean(axis=1)[:,np.newaxis])/rf_long.std(axis=1)[:,np.newaxis]
#



rf_short_anom = rf_short2 < -1.0
rf_long_anom = rf_long2 < -1.0
rf_short_fre = rf_short_anom.sum(axis=1)
rf_long_fre = rf_long_anom.sum(axis=1)



DataDict = {}
YM = [str(ym) for ym in range(2003,2019)]
for index,ym in enumerate(YM):
    DataDict["ShortEVI"+ym] = evi_short[:,index]
    DataDict["ShortRF" + ym] = rf_short[:, index]
    DataDict["ShortNSIF" + ym] = nsif_short[:, index]
    DataDict["ShortGSIF" + ym] = gsif_short[:, index]
    DataDict["ShortPVI" + ym] = spvi[:, index]

    DataDict["LongEVI" + ym] = evi_long[:, index]
    DataDict["LongRF" + ym] = rf_long[:, index]
    DataDict["LongNSIF" + ym] = nsif_long[:, index]
    DataDict["LongGSIF" + ym] = gsif_long[:, index]
    DataDict["LongPVI" + ym] = lpvi[:, index]


DataDict["ShortFre"] = rf_short_fre
DataDict["LongFre"] = rf_long_fre
df = pd.DataFrame(DataDict)
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask70_Mean.csv"
df.to_csv(outpath, index=False)

