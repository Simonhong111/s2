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
def monthAnalysis(path,monthtype,var1,year):
    if monthtype == "Short":
        months = [2,3,4,5]
    if monthtype == "Long":
        months =[6,7,8,9]
    data = pd.read_csv(path)
    key1 = [var1+str(year)+str(months[0]).zfill(2),var1+str(year)+str(months[1]).zfill(2),
            var1+str(year)+str(months[2]).zfill(2),var1+str(year)+str(months[3]).zfill(2)]


    Var1List = data[key1].to_numpy().mean(axis=1)

    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    return Var1List,PVIList



path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\CMIP\Cmip_PolyGonAgg_Mask50.csv"
Var1List,  PVIList = monthAnalysis(path, "Short", "RF", 2010)
row_num = Var1List.shape[0]
col_num = 13


rf_short = np.zeros(shape=(row_num,col_num),dtype=np.float)
spvi = np.zeros(shape=(row_num,col_num),dtype=np.float)

rf_long = np.zeros(shape=(row_num,col_num),dtype=np.float)
lpvi = np.zeros(shape=(row_num,col_num),dtype=np.float)

for id,Year in enumerate(range(2006,2019)):
    Var1List, PVIList = monthAnalysis(path, "Short", "RF", Year)

    rf_short[:,id] = Var1List
    spvi[:,id] = PVIList

    Var1List,  PVIList = monthAnalysis(path, "Long","RF", Year)

    rf_long[:, id] = Var1List
    lpvi[:, id] = PVIList



#
rf_short2 = (rf_short - rf_short.mean(axis=1)[:,np.newaxis])/rf_short.std(axis=1)[:,np.newaxis]
rf_long2 = (rf_long - rf_long.mean(axis=1)[:,np.newaxis])/rf_long.std(axis=1)[:,np.newaxis]
#



rf_short_anom = rf_short2 < -1.0
rf_long_anom = rf_long2 < -1.0
rf_short_fre = rf_short_anom.sum(axis=1)
rf_long_fre = rf_long_anom.sum(axis=1)



DataDict = {}
YM = [str(ym) for ym in range(2006,2019)]
for index,ym in enumerate(YM):

    DataDict["ShortRF" + ym] = rf_short[:, index]
    DataDict["ShortPVI" + ym] = spvi[:, index]


    DataDict["LongRF" + ym] = rf_long[:, index]
    DataDict["LongPVI" + ym] = lpvi[:, index]


DataDict["ShortFre"] = rf_short_fre
DataDict["LongFre"] = rf_long_fre
df = pd.DataFrame(DataDict)
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\CMIP\Cmip_PolyGonAgg_Mask50_Mean.csv"
df.to_csv(outpath, index=False)

