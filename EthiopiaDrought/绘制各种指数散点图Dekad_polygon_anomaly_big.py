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

def readVariableCSV(path,croptype,monthtype,year):
    if monthtype == "Short":
        months = [2,3,4,5]
    if monthtype == "Long":
        months = [6,7,8,9]

    cropkey = croptype+str(year)
    areakey = croptype[:-3]+"AREA"+str(year)
    rfkey = ["RF{}{}".format(str(year),str(m).zfill(2)) for m in months]
    pvikey = "{}PVI{}".format(monthtype,year)


    data = pd.read_csv(os.path.join(path, croptype + "_VegIndex_Big.csv"))

    Crop = data[cropkey].to_numpy()
    Area = data[areakey].to_numpy()
    print("*",Crop[0])
    print(Area[0])
    SCrop = Crop*Area
    print("***",SCrop[0],Crop[0]*Area[0])
    WCrop = SCrop.sum()/Area.sum()
    RF = data[rfkey].to_numpy().mean(axis=1).mean()
    PVI = data[pvikey].to_numpy().mean()



    return  WCrop,RF,PVI


def Anomly(path,croptype,monthtype):



    CropS = []
    RFS = []
    PVIS = []


    for i in range(7):
        Crop, RF, PVI = readVariableCSV(path, croptype, monthtype, i+2010)
        CropS.append(Crop)
        RFS.append(RF)
        PVIS.append(PVI)

    CropS = np.array(CropS)
    RFS = np.array(RFS)
    PVIS = np.array(PVIS)

    plt.scatter(CropS,RFS)
    plt.show()

    CropAnom = (CropS - CropS.mean())/CropS.std()
    RFAnom = (RFS - RFS.mean()) / RFS.std()
    PVIAnom = (PVIS - PVIS.mean()) / PVIS.std()


    slope, intercept, r_value, p_value, std_err = stats.linregress(RFAnom, PVIAnom)
    print("**",monthtype,croptype,r_value)


    return CropAnom,RFAnom,PVIAnom





DataType ="Polygon Mask"
path = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase"

CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH","CropAve"]


for croptype in CTp:
    fig = plt.figure(figsize=(15, 10))

    plt.title("2010-2016 Anomaly" + " For {}\n".format(croptype), fontsize=16)
    plt.xticks([])
    plt.yticks([])
    for m in range(1, 5):

        ax = fig.add_subplot(2, 2, m)

        if m == 1:
            titlename = "ShortRains RF Anom vs Yield Anom"
            print(m)
            CropAnom, RFAnom, PVIAnom = Anomly(path, croptype, "Short")
            Var1List, Var2List = RFAnom, CropAnom
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "RF Anomaly"


        if m == 2:
            titlename = "ShortRains PVI Anom vs Yield Anom"
            print(m)
            CropAnom, RFAnom, PVIAnom = Anomly(path, croptype, "Short")
            Var1List, Var2List = PVIAnom, CropAnom
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "PVI Anomaly"

        if m == 3:
            titlename = "LongRains RF Anom vs Yield Anom"
            print(m)
            CropAnom, RFAnom, PVIAnom = Anomly(path, croptype, "Long")
            Var1List, Var2List = RFAnom, CropAnom
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "RF Anomaly"

        if m == 4:
            titlename = "LongRains PVI Anom vs Yield Anom"
            print(m)
            CropAnom, RFAnom, PVIAnom = Anomly(path, croptype, "Long")
            Var1List, Var2List = PVIAnom, CropAnom
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "PVI Anomaly"

        scatter = ax.scatter(Var1List, Var2List)
        titles = titlename
        ax.set_title(titles, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Yield Anomaly")
        if m <= 2:
            ax.text(1, 0, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            ax.text(1, -0.2, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
        else:

            ax.text(-1.5, 0, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            ax.text(-1.5, -0.2, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
        for idx,_ in enumerate(Var1List):
            ax.text(Var1List[idx],Var2List[idx],str(idx+2010)[2:])

    fig.tight_layout()  # 调整整体空白
    path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase\Image10\Anom", croptype  + "Anom_AREA.jpg")
    plt.savefig(path2)
    plt.close()
    # plt.show()


# TestAnomly(path,"WHEATOPH","Short")