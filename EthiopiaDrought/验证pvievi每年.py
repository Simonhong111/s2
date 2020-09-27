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

def ExtractRFCSV(path,Vartype,Monthtype,defval=-9999):
    RF = []
    EVI = []

    CropType = os.path.basename(path).split("_")[0]
    if Monthtype == "short":
        Month = [2, 3, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))
    if Monthtype == "long":
        Month = [6, 7, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))


    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ctemp = []
            vartemp = []
            flag = False
            for ym in YM:
                if float(row[Vartype+ym]) == defval:
                    flag = True
            if flag:
                continue

            for year in range(2010,2017):
                ctemp.append(float(row[CropType+str(year)]))
            CropValue.append(ctemp)

            for ym in YM:
                vartemp.append(float(row[Vartype+ym]))
            VarValue.append(vartemp)
    number = len(VarValue)
    CropValue = np.array(CropValue)

    VarValue = np.array(VarValue).reshape((number,7,4))
    VarValueNew = VarValue.mean(axis=2)
    VarValueMean = VarValueNew.mean(axis=1)
    VarValueStd = VarValueNew.std(axis=1)
    CropMean = CropValue.mean(axis=1)
    CropStd = CropValue.std(axis=1)

    for id in range(number):
        CropValue[id,:] = (CropValue[id,:] - CropMean[id])/CropStd[id]
        VarValueNew[id,:] = (VarValueNew[id,:] - VarValueMean[id])/VarValueStd[id]

    drought = np.zeros(shape=(number,1),dtype=np.float)

    for id in range(number):
        drought[id] = np.where(VarValueNew[id,:] <= -1.0)[0].shape[0]

    return CropValue,VarValueNew,drought

def ExtractEVICSV(path,Vartype,Monthtype,defval=-3000):
    CropValue = []

    VarValue = []
    RFValue = []
    CropType = os.path.basename(path).split("_")[0]
    if Monthtype == "short":
        Month = [2, 3, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))
    if Monthtype == "long":
        Month = [6, 7, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))


    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ctemp = []
            vartemp = []
            rftemp = []
            flag = False
            for ym in YM:
                if float(row[Vartype+ym]) == defval or float(row["RF"+ym]) == -9999:
                    flag = True
            if flag:
                continue

            for year in range(2010,2017):
                ctemp.append(float(row[CropType+str(year)]))
            CropValue.append(ctemp)


            for ym in YM:
                vartemp.append(float(row[Vartype+ym]))
                rftemp.append(float(row["RF"+ym]))
            VarValue.append(vartemp)
            RFValue.append(rftemp)
    number = len(VarValue)
    CropValue = np.array(CropValue)
    CropMean = CropValue.mean(axis=1)
    CropStd = CropValue.std(axis=1)

    VarValue = np.array(VarValue).reshape((number,7,4))
    VarValueNew = VarValue.mean(axis=2)
    VarValueMean = VarValueNew.mean(axis=1)
    VarValueStd = VarValueNew.std(axis=1)

    RFValue = np.array(RFValue).reshape((number,7,4))
    RFValueNew = RFValue.mean(axis=2)
    RFValueMean = RFValueNew.mean(axis=1)
    RFValueStd = RFValueNew.std(axis=1)


    for id in range(number):
        CropValue[id,:] = (CropValue[id,:] - CropMean[id])/CropStd[id]
        VarValueNew[id,:] = (VarValueNew[id,:] - VarValueMean[id])/VarValueStd[id]
        RFValueNew[id,:] = (RFValueNew[id,:] - RFValueMean[id])/RFValueStd[id]
    drought = np.zeros(shape=(number,1),dtype=np.float)

    for id in range(number):
        drought[id] = np.where(RFValueNew[id,:] <= -1.0)[0].shape[0]
    return CropValue,VarValueNew,drought


def ExtractPVICSV(path,Vartype,Monthtype,defval=-9999):
    CropValue = []

    VarValue = []
    RFValue = []
    CropType = os.path.basename(path).split("_")[0]
    if Monthtype == "short":
        Month = [2, 3, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))
    if Monthtype == "long":
        Month = [6, 7, 4, 5]
        Year = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
        YM = []
        for year in Year:
            for m in Month:
                YM.append(str(year) + str(m).zfill(2))


    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ctemp = []
            vartemp = []
            rftemp = []
            flag = False
            for y in range(2010,2017):
                if float(row[Vartype + str(y)]) == defval:
                    flag = True
            for ym in YM:
                if float(row["RF"+ym]) == -9999:
                    flag = True
            if flag:
                continue

            for year in range(2010,2017):
                ctemp.append(float(row[CropType+str(year)]))
                vartemp.append(float(row[Vartype+str(year)]))
            CropValue.append(ctemp)
            VarValue.append(vartemp)


            for ym in YM:

                rftemp.append(float(row["RF"+ym]))

            RFValue.append(rftemp)
    number = len(VarValue)
    CropValue = np.array(CropValue)
    CropMean = CropValue.mean(axis=1)
    CropStd = CropValue.std(axis=1)

    VarValue = np.array(VarValue)
    VarValueMean = VarValue.mean(axis=1)
    VarValueStd = VarValue.std(axis=1)

    RFValue = np.array(RFValue).reshape((number,7,4))
    RFValueNew = RFValue.mean(axis=2)
    RFValueMean = RFValueNew.mean(axis=1)
    RFValueStd = RFValueNew.std(axis=1)


    for id in range(number):
        CropValue[id,:] = (CropValue[id,:] - CropMean[id])/CropStd[id]
        VarValue[id,:] = (VarValue[id,:] - VarValueMean[id])/VarValueStd[id]
        RFValueNew[id,:] = (RFValueNew[id,:] - RFValueMean[id])/RFValueStd[id]
    drought = np.zeros(shape=(number,1),dtype=np.float)

    for id in range(number):
        drought[id] = np.where(RFValueNew[id,:] <= -1.0)[0].shape[0]
    return CropValue,VarValue,drought


path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\BARLEYOPH_From2010-2016_V2.csv"
YM = []
for year in range(2010,2017):
    for month in range(1,13):
        YM.append(str(year)+str(month).zfill(2))

R = []
for ym in YM:

    VarList1,VarList2 = ExtractRFCSV(path,"RF","EVI",ym)
    slope, intercept, r_value, p_value, std_err = stats.linregress(VarList1,VarList2)
    title = "Time:{} Correlation R between Rainfall and EVI {}".format(ym,r_value)
    print(title)
    R.append(r_value)

R2=[]
for ym in YM:

    VarList1,VarList2 = ExtractPVICSV(path,"EVI","PVI",ym)
    slope, intercept, r_value, p_value, std_err = stats.linregress(VarList1,VarList2)
    title = "Time:{} Correlation R between PVI and EVI {}".format(ym,r_value)
    print(title)
    R2.append((-1)*r_value)


print(len(R),len(R2))

plt.bar(range(len(R)),R,label="evi_rainfall")
plt.plot(range(len(R2)),R2,"r*",label='evi_pvi')
plt.plot(range(len(R2)),R2,"r")
plt.xticks(range(len(R)),YM)
plt.ylabel("R value")
plt.xlabel("time year&month")

for tick in plt.gca().get_xticklabels():

    tick.set_rotation(90)
plt.legend()
plt.show()


