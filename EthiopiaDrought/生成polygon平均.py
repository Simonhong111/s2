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

def readCropCSV():
    path = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"

    CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]

    keys = []
    wkeys = []
    for year in range(2010,2017):
        for ctp in CTp:
            keys.append(ctp+str(year))
            wkeys.append(ctp[:-3]+"AREA"+str(year))

    print(wkeys)

    data = pd.read_csv(path)
    ID = data["_ID"].to_numpy()
    crop_subset = data[keys].to_numpy()
    area_subset = data[wkeys].to_numpy()

    invalidRow = []
    for i in range(7):
        subset = crop_subset[:,(i*5):(i*5+5)]

        for row in range(subset.shape[0]):
            NAMask = np.isnan(subset[row])*1
            if NAMask.sum() == 0:
                invalidRow.append(row)
    invalidRow = list(set(invalidRow))
    validRow = []
    for i in  range(crop_subset.shape[0]):
        if i not in invalidRow:
            validRow.append(i)

    crop_subset = crop_subset[validRow,:]
    area_subset = area_subset[validRow,:]
    ID = ID[validRow]

    cropmask = np.isnan(crop_subset)
    crop_subset[cropmask] = 0
    area_subset[cropmask] = 0
    areamask = np.isnan(area_subset)
    area_subset[areamask] = 0
    crop_subset[areamask] = 0


    DataDict = {}

    for year in range(2010, 2017):
        col_st = year -2010

        temp = crop_subset[:,(col_st*5):(col_st*5+5)]*area_subset[:,(col_st*5):(col_st*5+5)]
        DataDict["CropAve" + str(year)] = temp.sum(axis=1)/area_subset[:,(col_st*5):(col_st*5+5)].sum(axis=1)
        for idx,ctp in enumerate(CTp):
            DataDict[ctp + str(year)] = crop_subset[:,(col_st*5):(col_st*5+5)][:,idx]
            DataDict[ctp[:-3] +"AREA"+ str(year)] = area_subset[:, (col_st * 5):(col_st * 5 + 5)][:, idx]


    DataDict["ID"] = ID



    df = pd.DataFrame(DataDict)
    outpath = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase\{}_VegIndex.csv".format("CropAve")
    df.to_csv(outpath, index=False)


# readCropCSV()

DataType ="Crop Mask"
path = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase"


# 验证
def Test():
    path = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"

    CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]

    keys = []
    wkeys = []
    for year in range(2010,2017):
        for ctp in CTp:
            keys.append(ctp+str(year))
            wkeys.append(ctp[:-3]+"AREA"+str(year))


    data = pd.read_csv(path)
    ID = data["_ID"].to_numpy()
    crop_subset = data[keys].to_numpy()
    area_subset = data[wkeys].to_numpy()

    selecRow = np.where(ID == 40930)

    print(crop_subset[selecRow[0],5:10])
    print(area_subset[selecRow[0], 5:10])

Test()

a = np.array([        0, 25.39999962 ,      0, 28.00000381, 12.69999313])
b = np.array([ 0.      ,   13.80362201,  0.,          0.754556  ,  2.37596498])
c = a*b
print(c.sum()/b.sum())


