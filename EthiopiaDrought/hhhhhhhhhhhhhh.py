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
    CropValue = []
    VarValue = []
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
def ExtractNEWSIFCSV(path,Vartype,Monthtype,defval=-9999):
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
def ExtractGOSIFCSV(path,Vartype,Monthtype,defval=32766):
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

                if float(row[Vartype+ym]) >= defval or float(row["RF"+ym]) == -9999:
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

def ExtractETCSV(path,Vartype,Monthtype,defval=-999):
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
def ExtractSMCSV(path,Vartype,Monthtype,defval=-2):
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
directory = r"D:\Cornell\EthiopianDrought\CropCSV\Crop"
# croptype =  ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
# path = os.path.join(directory,croptype[0]+"_From2010-2016.csv")
# CropValue,VarValue,drought = ExtractGOSIFCSV(path,"GOSIF","long",32766)
#





CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
# cropclass = CTp[4]
# mMonthType = "long"
# year = 6
Years =["2010","2011","2012","2013","2014","2015","2016"]

colors =['b','g','orange','brown','r']
for cropclass in CTp:

    path = os.path.join(directory, cropclass + "_From2010-2016_V2.csv")
    for mMonthType in ["short","long"]:

        for year in range(7):
            fig = plt.figure(figsize=(15, 10))
            plt.title(cropclass + " " + Years[year] + " " + mMonthType + " rains" + "\n",fontsize=16)
            plt.ylabel("Yield Anomaly",fontsize=16)
            for m in range(1, 8):

                ax = fig.add_subplot(3, 3, m)

                if m == 1:

                    varname = "RainFall r:"
                    xlabel = "RainFall Anomaly"
                    print(varname)
                    CropValue, VarValue, drought = ExtractRFCSV(path, "RF", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                if m == 2:
                    varname = "EVI r:"
                    print(varname)
                    xlabel = "EVI Anomaly"
                    CropValue, VarValue, drought = ExtractEVICSV(path, "EVI", mMonthType)
                    print(CropValue[0][0])
                    droughtID = list(set(list(drought[:, 0])))
                if m == 3:
                    varname = "NEWSIF r:"
                    print(varname)
                    xlabel = "NEWSIF Anomaly"
                    CropValue, VarValue, drought = ExtractNEWSIFCSV(path, "NEWSIF", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                if m == 4:
                    varname = "GOSIF r:"
                    print(varname)
                    xlabel = "GOSIF Anomaly"
                    CropValue, VarValue, drought = ExtractGOSIFCSV(path, "GOSIF", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                if m == 5:
                    varname = "ET r:"
                    print(varname)
                    xlabel = "ET Anomaly"
                    CropValue, VarValue, drought = ExtractETCSV(path, "ET", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                if m == 6:
                    varname = "SM r:"
                    print(varname)
                    xlabel = "SM Anomaly"
                    CropValue, VarValue, drought = ExtractSMCSV(path, "SM", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                if m == 7:
                    varname = "PVI r:"
                    print(varname)
                    xlabel = "PVI Anomaly"
                    CropValue, VarValue, drought = ExtractPVICSV(path, "PVI", mMonthType)
                    droughtID = list(set(list(drought[:, 0])))
                mcolor = []
                for i in range(drought.shape[0]):
                    if drought[i, 0] == 0:
                        mcolor.append(colors[0])
                    if drought[i, 0] == 1:
                        mcolor.append(colors[1])
                    if drought[i, 0] == 2:
                        mcolor.append(colors[2])
                    if drought[i, 0] == 3:
                        mcolor.append(colors[3])
                    if drought[i, 0] == 4:
                        mcolor.append(colors[4])
                    if drought[i, 0] == 5:
                        mcolor.append(colors[5])
                colors = np.array(colors)

                scatter = ax.scatter(VarValue[:, year], CropValue[:, year], c=mcolor)
                titles = "R for freq "
                for freq in droughtID:
                    titles = titles+str(int(freq))+" "
                titles=titles+":"
                for freq in droughtID:
                    mask = np.where(drought[:, 0] == freq)
                    sample = mask[0].shape[0]
                    slope, intercept, r_value, p_value, std_err = stats.linregress(VarValue[:, year][mask],
                                                                                   CropValue[:, year][mask])
                    titles = titles + str(round(r_value, 2)) + ","
                ax.tick_params(labelsize=8)
                ax.set_title(titles, fontsize=14)
                ax.set_xlabel(xlabel,fontsize=12)
                # if m == 7:
                #     legend1 = ax.legend(*scatter.legend_elements(),
                #             loc="upper left", bbox_to_anchor=(1.0,1.0),title="Drought Frequency")
                #     ax.add_artist(legend1)
            ax = fig.add_subplot(3, 3, 8)
            ax.scatter([1, 1, 1, 1, 1], [2.5, 2, 1.5, 1, 0.5], c=['b', 'g', 'orange', 'brown', 'r'])
            ax.text(1.005, 2.485,"Frequency 0")
            ax.text(1.005, 1.985, "Frequency 1")
            ax.text(1.005, 1.485, "Frequency 2")
            ax.text(1.005, 0.985, "Frequency 3")
            ax.text(1.005, 0.485, "Frequency 4")
            fig.tight_layout()  # 调整整体空白
            path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\image",
                                 cropclass + mMonthType + Years[year] + ".jpg")
            print(path2)
            plt.savefig(path2)
            plt.close(fig)

# plt.show()
