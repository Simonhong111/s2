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


def ExtractRFCSV(path,Var1,Var2,YM,defv1=-9999,defv2=-3000):
    VarList1 = []
    VarList2 = []

    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row[Var1+YM] == defv1 or row[Var2+YM] == defv2:
                continue
            VarList1.append(float(row[Var1+YM]))
            VarList2.append(float(row[Var2+YM]))

    return VarList1,VarList2
def ExtractPVICSV(path,Var1,Var2,YM,defv1=-3000,defv2=-9999):
    VarList1 = []
    VarList2 = []

    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row[Var1+YM] == defv1 or row[Var2+YM[0:4]] == defv2:
                continue
            VarList1.append(float(row[Var1+YM]))
            VarList2.append(float(row[Var2+YM[0:4]]))

    return VarList1,VarList2


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


