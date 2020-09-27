import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\10day\PolyGonAgg_Mask70.csv"

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

    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    return Var1List,Var2List,PVIList
# monthAnalysis(path,"Short","EVI","RF","2010")

DataType ="Crop Mask"
Years = [str(year) for year in range(2003,2019)]
for Year in Years:
    fig = plt.figure(figsize=(15, 10))
    mMonthType = "Long"
    plt.title(Year + " For {}\n".format(DataType), fontsize=16)
    plt.xticks([])
    plt.yticks([])
    for m in range(1, 7):

        ax = fig.add_subplot(2, 3, m)
        if m == 1:
            titlename = "{} ShortRains EVI vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Short", "EVI", "RF", Year)
            xlabel = "EVI"

            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        if m == 2:
            titlename = "{} ShortRains GOSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Short", "GSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "GOSIF"
        if m == 3:
            titlename = "{} ShortRains NewSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Short", "NSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "NEWSIF"
        if m == 4:
            titlename = "{} LongRains EVI vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Long", "EVI", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "EVI"
        if m == 5:
            titlename = "{} LongRains GOSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Long", "GSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "GOSIF"
        if m == 6:
            titlename = "{} LongRains NSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List = monthAnalysis(path, "Long", "NSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
            xlabel = "NEWSIF"

        scatter = ax.scatter(Var1List, Var2List)
        titles = titlename
        ax.set_title(titles, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("RainFall")
        if m <= 3:
            ax.text(0.1, 190, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            # ax.text(0.1, 180, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
        else:

            ax.text(0.5, 380, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            ax.text(0.5, 360, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")

    fig.tight_layout()  # 调整整体空白
    path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\10day\Image10", DataType + Year + "RF70.jpg")
    plt.savefig(path2)
    plt.close()
    # plt.show()

#
#
