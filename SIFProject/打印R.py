import os
import glob
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv
from scipy import stats
from matplotlib import font_manager as fm, rcParams
import pandas as pd
global Va1,Va2,Va3,Va4,Va5,Va6,Va7,Va8

Va1,Va2,Va3,Va4,Va5,Va6,Va7,Va8,V0 = [],[],[],[],[],[],[],[],[]

fpath = os.path.join(rcParams["datapath"], r"C:\Users\Administrator\Desktop\finalexcel_hzm\times.ttf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(20)
prop.set_weight('bold')
proplabel = prop.copy()
proplabel.set_size(18)
proptext = prop.copy()
proptext.set_size(16)


# slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
# mpl.rcParams['font.sans-serif'] = ['FangSong']
# mpl.rcParams['axes.unicode_minus'] = False
spefied_year = 2009

def reg(path, V1, V2, key):
    Data1 = []
    Data2 = []
    Date = []
    YDate = []
    with open(path, encoding='UTF-8') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            Data1.append(float(row[V1]))
            Data2.append(float(row[V2]))
            Date.append(float(row["month"]))
            YDate.append(float(row["year"]))
    Date = np.array(Date)
    YDate = np.array(YDate)
    Data1 = np.array(Data1)
    Data2 = np.array(Data2)

    if key == 0:
        #   非汛期
        # date 月 Ydate 年
        extent = (((Date < 7) & (Date > 0)) | (Date > 8)) & (YDate >= spefied_year)
    if key == 1:
        print(key)
        # 汛前
        extent = (Date > 3) & (Date < 7) & (YDate >= spefied_year)
    if key == 2:
        # 汛期
        extent = ((Date > 6) & (Date <= 9)) & (YDate >= spefied_year)
    if key == 3:
        # 蓄水期
        extent = ((Date > 9) & (Date <= 12)) & (YDate >= spefied_year)

    if key == 4:
        # 枯水期
        extent = (((Date > 0) & (Date <= 6)) | (Date == 12)) & (YDate >= spefied_year)
    if key == 5:
        # 枯水期
        extent = (Date > 0) & (YDate >= spefied_year)


    w = np.where(extent)

    slope, intercept, r_value, p_value, std_err = stats.linregress(Data1[w], Data2[w])

    if str(round(p_value, 3)).__len__() < 5:

        p_value = '0.000'
    else:
        p_value = str(round(p_value, 3))
    R2 = round(r_value**2,3)
    intercept = round(intercept,3)
    slope = round(slope,6)

    return slope, intercept, R2, p_value, Data1[w], Data2[w]




def prt(dir, V1,V2,key=5):
    rsts = glob.glob(os.path.join(dir, "*10days.csv"))
    for rst in rsts:
        # V1 = "water"
        # V2 = "avg"
        slope, intercept, r_squa, p_value, Data1, Data2 = reg(rst, V1, V2, key=key)
        if "10" in os.path.basename(rst)[5:]:
            interval = 10
        if "30" in os.path.basename(rst)[5:]:
            interval = 30

        if key == 6:
            #   非汛期
            mask = "非汛期"
        if key == 7:
            #   非汛期
            mask = "汛期与蓄水"
        if key == 1:
            # 汛前
            mask = "汛前"
        if key == 2:
            # 汛期
            mask = "汛期"
        if key == 3:
            # 蓄水期
            mask = "蓄水期"

        if key == 4:
            # 枯水期
            mask = "消落期"
        if key == 5:
            # 枯水期
            mask = "全年"

        V11 = ""
        if V1 == "water":
            V11 = "水位"
        if V1 == "preciption":
            V11 = "降雨"
        if V1 == "avg":
            V11 = "温度"

        V22 = ""
        if V2 == "preciption":
            V22 = u"降雨"
        if V2 == "avg":
            V22 = u"温度"
        if V2 =="ndvi":
            V22 = u"植被"

        V0.append(os.path.basename(rst).split("_")[1])
        Va1.append(V11)
        Va2.append(V22)
        Va3.append(mask)
        Va4.append(interval)
        Va5.append(slope)
        Va6.append(intercept)
        Va7.append(r_squa)
        Va8.append(p_value)



Name1 = ["water","water","preciption","avg"]
Name2 = ["preciption","avg","ndvi","ndvi"]

for idx,item in enumerate(Name1):
    #prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2",item,Name2[idx], 1)
    prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item,Name2[idx],2)
    prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item,Name2[idx],3)
    prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item,Name2[idx],4)
    prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item,Name2[idx],5)
    #prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item, Name2[idx], 6)
    #prt(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Result2", item, Name2[idx], 7)

dataframe = pd.DataFrame({"Distance":V0,'Variable1': Va1, 'Variable2': Va2, "stage": Va3, "interval": Va4, "slope": Va5, "intercept": Va6,"r2": Va7, "p-value": Va8})

    # 将DataFrame存储为csv,index表示是否显示行名，default=True
dataframe.to_csv(r"C:\Users\Administrator\Desktop\finalexcel_hzm\10_dongying_2009_2015.csv", index=False, sep=',',encoding="gbk")







