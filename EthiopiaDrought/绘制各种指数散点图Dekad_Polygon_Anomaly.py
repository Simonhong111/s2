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
    rfkey = ["RF{}{}".format(str(year),str(m).zfill(2)) for m in months]
    pvikey = "{}PVI{}".format(monthtype,year)
    evikey = ["EVI{}{}".format(str(year), str(m).zfill(2)) for m in months]
    gsifkey = ["GOSIF{}{}".format(str(year), str(m).zfill(2)) for m in months]
    nsifkey = ["NEWSIF{}{}".format(str(year), str(m).zfill(2)) for m in months]

    data = pd.read_csv(os.path.join(path, croptype + "_VegIndex.csv"))

    Crop = data[cropkey].to_numpy()
    RF = data[rfkey].to_numpy().mean(axis=1)
    PVI = data[pvikey].to_numpy()
    EVI = data[evikey].to_numpy().mean(axis=1)
    GSIF = data[gsifkey].to_numpy().mean(axis=1)
    NSIF = data[nsifkey].to_numpy().mean(axis=1)


    return  Crop,RF,PVI,EVI,GSIF,NSIF


def Anomly(path,croptype,monthtype):

    Crop,RF,PVI,EVI,GSIF,NSIF = readVariableCSV(path, croptype, monthtype, 2010)
    Num = Crop.shape[0]
    CropS = np.zeros((Num,7),dtype=np.float)
    RFS = np.zeros((Num, 7), dtype=np.float)
    PVIS = np.zeros((Num, 7), dtype=np.float)
    EVIS = np.zeros((Num, 7), dtype=np.float)
    GSIFS = np.zeros((Num, 7), dtype=np.float)
    NSIFS = np.zeros((Num, 7), dtype=np.float)

    for i in range(7):
        Crop, RF, PVI, EVI, GSIF, NSIF = readVariableCSV(path, croptype, monthtype, i+2010)
        CropS[:,i] = Crop
        RFS[:,i] = RF
        PVIS[:,i] = PVI
        EVIS[:,i] = EVI
        GSIFS[:,i] = GSIF
        NSIFS[:,i] = NSIF

    CropAnom = (CropS - CropS.mean(axis=0))/CropS.std(axis=0)
    RFAnom = (RFS - RFS.mean(axis=0)) / RFS.std(axis=0)
    PVIAnom = (PVIS - PVIS.mean(axis=0)) / PVIS.std(axis=0)
    EVIAnom = (EVIS - EVIS.mean(axis=0)) / EVIS.std(axis=0)
    GSIFAnom = (GSIFS - GSIFS.mean(axis=0)) / GSIFS.std(axis=0)
    NSIFAnom = (NSIFS - NSIFS.mean(axis=0)) / NSIFS.std(axis=0)
    return CropAnom,RFAnom,PVIAnom,EVIAnom,GSIFAnom,NSIFAnom





DataType ="Polygon Mask"
path = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase"

CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH","CropAve"]
#
# Years = [year for year in range(2010,2017)]
# for croptype in CTp:
#     for Year in Years:
#         fig = plt.figure(figsize=(15, 10))
#
#         plt.title(str(Year) + " For {}\n".format(DataType), fontsize=16)
#         plt.xticks([])
#         plt.yticks([])
#         for m in range(1, 9):
#
#             ax = fig.add_subplot(2, 4, m)
#             if m == 1:
#                 titlename = "{} ShortRains EVI Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Short")
#                 Var1List,Var2List = EVIAnom[:,Year-2010],CropAnom[:,Year-2010]
#                 xlabel = "EVI"
#
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#             if m == 2:
#                 titlename = "{} ShortRains RF Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Short")
#                 Var1List, Var2List = RFAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "RF"
#             if m == 3:
#                 titlename = "{} ShortRains NewSIF Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Short")
#                 Var1List, Var2List = NSIFAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "NEWSIF"
#             if m == 4:
#                 titlename = "{} ShortRains PVI Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Short")
#                 Var1List, Var2List = PVIAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "PVI"
#             if m == 5:
#                 titlename = "{} LongRains EVI Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Long")
#                 Var1List, Var2List = EVIAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "EVI"
#             if m == 6:
#                 titlename = "{} LongRains RF Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Long")
#                 Var1List, Var2List = RFAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "RF"
#             if m == 7:
#                 titlename = "{} LongRains NSIF Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Long")
#                 Var1List, Var2List = NSIFAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "NEWSIF"
#             if m == 8:
#                 titlename = "{} LongRains PVI Anom vs Yield Anom".format(Year)
#                 CropAnom, RFAnom, PVIAnom, EVIAnom, GSIFAnom, NSIFAnom = Anomly(path, croptype, "Long")
#                 Var1List, Var2List = PVIAnom[:, Year - 2010], CropAnom[:, Year - 2010]
#                 slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#                 xlabel = "PVI"
#
#             scatter = ax.scatter(Var1List, Var2List)
#             titles = titlename
#             ax.set_title(titles, fontsize=10)
#             ax.set_xlabel(xlabel)
#             ax.set_ylabel("Yield")
#             if m <= 4:
#                 ax.text(0.1, 0.3, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#                 # ax.text(0.1, 180, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#             else:
#
#                 ax.text(0.1, 0.1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#                 # ax.text(0.1, 0.1, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#
#         fig.tight_layout()  # 调整整体空白
#         path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase\Image10\Anom", croptype + str(Year) + "Anom.jpg")
#         plt.savefig(path2)
#         plt.close()
#         # plt.show()

#
#

def TestAnomly(path,croptype,monthtype):

    Crop,RF,PVI,EVI,GSIF,NSIF = readVariableCSV(path, croptype, monthtype, 2010)
    Num = Crop.shape[0]
    CropS = np.zeros((Num,7),dtype=np.float)
    RFS = np.zeros((Num, 7), dtype=np.float)
    PVIS = np.zeros((Num, 7), dtype=np.float)
    EVIS = np.zeros((Num, 7), dtype=np.float)
    GSIFS = np.zeros((Num, 7), dtype=np.float)
    NSIFS = np.zeros((Num, 7), dtype=np.float)

    for i in range(7):
        Crop, RF, PVI, EVI, GSIF, NSIF = readVariableCSV(path, croptype, monthtype, i+2010)
        CropS[:,i] = Crop
        RFS[:,i] = RF
        PVIS[:,i] = PVI
        EVIS[:,i] = EVI
        GSIFS[:,i] = GSIF
        NSIFS[:,i] = NSIF

    CropAnom = (CropS - CropS.mean(axis=0))/CropS.std(axis=0)
    print(CropS[2,:],CropS.mean(axis=0)[2])
    print((CropS[2, :]- CropS.mean(axis=0))/CropS.std(axis=0))
    print((CropS - CropS.mean(axis=0))[2])
    print(CropAnom[2])

    RFAnom = (RFS - RFS.mean(axis=0)) / RFS.std(axis=0)
    PVIAnom = (PVIS - PVIS.mean(axis=0)) / PVIS.std(axis=0)
    EVIAnom = (EVIS - EVIS.mean(axis=0)) / EVIS.std(axis=0)
    GSIFAnom = (GSIFS - GSIFS.mean(axis=0)) / GSIFS.std(axis=0)
    NSIFAnom = (NSIFS - NSIFS.mean(axis=0)) / NSIFS.std(axis=0)


# TestAnomly(path,"WHEATOPH","Short")