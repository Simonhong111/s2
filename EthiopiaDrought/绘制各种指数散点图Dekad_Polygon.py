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


DataType ="Polygon Mask"
path = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase"

CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH","CropAve"]
Years = [str(year) for year in range(2010,2017)]
for croptype in CTp:
    for Year in Years:
        fig = plt.figure(figsize=(15, 10))

        plt.title(Year + " For {}\n".format(DataType), fontsize=16)
        plt.xticks([])
        plt.yticks([])
        for m in range(1, 9):

            ax = fig.add_subplot(2, 4, m)
            if m == 1:
                titlename = "{} ShortRains EVI vs Yield".format(Year)
                Var2List,RF,PVI,Var1List,GSIF,NSIF = readVariableCSV(path,croptype,"Short",Year)
                xlabel = "EVI"

                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            if m == 2:
                titlename = "{} ShortRains RF vs Yield".format(Year)
                Var2List,Var1List,PVI,EVI,GSIF,NSIF = readVariableCSV(path,croptype,"Short",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "RF"
            if m == 3:
                titlename = "{} ShortRains NewSIF vs Yield".format(Year)
                Var2List,RF,PVI,EVI,GSIF,Var1List = readVariableCSV(path,croptype,"Short",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "NEWSIF"
            if m == 4:
                titlename = "{} ShortRains PVI vs Yield".format(Year)
                Var2List,RF,Var1List,EVI,GSIF,NSIF = readVariableCSV(path,croptype,"Short",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "PVI"
            if m == 5:
                titlename = "{} LongRains EVI vs Yield".format(Year)
                Var2List,RF,PVI,Var1List,GSIF,NSIF = readVariableCSV(path,croptype,"Long",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "EVI"
            if m == 6:
                titlename = "{} LongRains RF vs Yield".format(Year)
                Var2List,Var1List,PVI,EVI,GSIF,NSIF = readVariableCSV(path,croptype,"Long",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "RF"
            if m == 7:
                titlename = "{} LongRains NSIF vs Yield".format(Year)
                Var2List,RF,PVI,EVI,GSIF,Var1List = readVariableCSV(path,croptype,"Long",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "NEWSIF"
            if m == 8:
                titlename = "{} LongRains PVI vs Yield".format(Year)
                Var2List,RF,Var1List,EVI,GSIF,NSIF = readVariableCSV(path,croptype,"Long",Year)
                slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
                xlabel = "PVI"

            scatter = ax.scatter(Var1List, Var2List)
            titles = titlename
            ax.set_title(titles, fontsize=10)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Yield")
            if m <= 4:
                ax.text(0.1, 0.3, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
                # ax.text(0.1, 180, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
            else:

                ax.text(0.1, 0.1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
                # ax.text(0.1, 0.1, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")

        fig.tight_layout()  # 调整整体空白
        path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase\Image10\Ave", croptype + Year + ".jpg")
        plt.savefig(path2)
        plt.close()
        # plt.show()

#
#

