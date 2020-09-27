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
def monthAnalysis(path):

    data = pd.read_csv(path)

    short_evi_name = ["Short"+"EVI"+str(year) for year in range(2007,2019)]
    long_evi_name = ["Long" + "EVI" + str(year) for year in range(2007, 2019)]
    short_rf_name = ["Short" + "RF" + str(year) for year in range(2007, 2019)]
    long_rf_name = ["Long" + "RF" + str(year) for year in range(2007, 2019)]
    short_pvi_name = ["Short" + "PVI" + str(year) for year in range(2007, 2019)]
    long_pvi_name = ["Long" + "PVI" + str(year) for year in range(2007, 2019)]
    short_gsif_name = ["Short" + "GSIF" + str(year) for year in range(2007, 2019)]
    long_gsif_name = ["Long" + "GSIF" + str(year) for year in range(2007, 2019)]
    short_nsif_name = ["Short" + "NSIF" + str(year) for year in range(2007, 2019)]
    long_nsif_name = ["Long" + "NSIF" + str(year) for year in range(2007, 2019)]
    short_Fre_name = ["Short" + "Fre"]
    long_Fre_name = ["Long" + "Fre" ]

    short_evi = data[short_evi_name].to_numpy()
    long_evi =data[long_evi_name].to_numpy()
    short_rf =data[short_rf_name].to_numpy()
    long_rf =data[long_rf_name].to_numpy()
    short_pvi =data[short_pvi_name].to_numpy()
    long_pvi =data[long_pvi_name].to_numpy()
    short_gsif =data[short_gsif_name].to_numpy()
    long_gsif =data[long_gsif_name].to_numpy()
    short_nsif =data[short_nsif_name].to_numpy()
    long_nsif =data[long_nsif_name].to_numpy()
    short_Fre =data[short_Fre_name].to_numpy()
    long_Fre = data[long_Fre_name].to_numpy()


    return short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif ,short_Fre ,long_Fre

path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\new\PolyGonAgg_Mask70_Mean.csv"

short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif ,short_Fre ,long_Fre=monthAnalysis(path)

Years = range(2007,2019)
Ylabel = [str(year)[-2:] for year in range(2007,2019)]
ShortFreSet = set(short_Fre.flatten().tolist())
LongFreSet = set(long_Fre.flatten().tolist())
print(ShortFreSet)
print(LongFreSet)

short_evi_anom = (short_evi.mean(axis=0) - short_evi.mean(axis=0).mean())/short_evi.mean(axis=0).std()
long_evi_anom = (long_evi.mean(axis=0) - long_evi.mean(axis=0).mean())/long_evi.mean(axis=0).std()
short_pvi_anom = (short_pvi.mean(axis=0) - short_pvi.mean(axis=0).mean())/short_pvi.mean(axis=0).std()
long_pvi_anom = (long_pvi.mean(axis=0) - long_pvi.mean(axis=0).mean())/long_pvi.mean(axis=0).std()
short_gsif_anom = (short_gsif.mean(axis=0) - short_gsif.mean(axis=0).mean())/short_gsif.mean(axis=0).std()
long_gsif_anom = (long_gsif.mean(axis=0) - long_gsif.mean(axis=0).mean())/long_gsif.mean(axis=0).std()
short_nsif_anom = (short_nsif.mean(axis=0) - short_nsif.mean(axis=0).mean())/short_nsif.mean(axis=0).std()
long_nsif_anom = (long_nsif.mean(axis=0) - long_nsif.mean(axis=0).mean())/long_nsif.mean(axis=0).std()
short_rf_anom = (short_rf.mean(axis=0) - short_rf.mean(axis=0).mean())/short_rf.mean(axis=0).std()
long_rf_anom = (long_rf.mean(axis=0) - long_rf.mean(axis=0).mean())/long_rf.mean(axis=0).std()
print(short_evi.mean(axis=0),short_evi.mean(axis=0).std())

fig = plt.figure(figsize=(15, 10))
plt.title("area-average anomaly of short rains, 2003 - 2018 \n")
plt.xticks([])
plt.yticks([])

ax1 = fig.add_subplot(3, 1, 1)
ax1.set_title("evi average anomaly")
ax1.grid(b=True)
ax1.plot(Years,short_evi_anom,label='evi')
ax1.plot(Years,short_pvi_anom,'r--',label='pvi')
ax1.plot(Years,short_rf_anom,'k--',label='rainfall')
ax1.set_xticks(Years)
ax1.set_xticklabels(Ylabel)
ax1.legend()

ax2 = fig.add_subplot(3, 1, 2)
ax2.set_title("\n\ngsif average anomaly")
ax2.grid(b=True)
ax2.plot(Years,short_gsif_anom,label='gosif')
ax2.plot(Years,short_pvi_anom,'r--',label='pvi')
ax2.plot(Years,short_rf_anom,'k--',label='rainfall')
ax2.set_xticks(Years)
ax2.set_xticklabels(Ylabel)
ax2.legend()


ax3 = fig.add_subplot(3, 1, 3)
ax3.set_title("\n\nnsif average anomaly")
ax3.grid(b=True)
ax3.plot(Years,short_nsif_anom,label='newsif')
ax3.plot(Years,short_pvi_anom,'r--',label='pvi')
ax3.plot(Years,short_rf_anom,'k--',label='rainfall')
ax3.set_xticks(Years)
ax3.set_xticklabels(Ylabel)
ax3.legend()

# fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0, hspace =0.3)#调整子图间距
fig = plt.figure(figsize=(15, 10))
plt.title("area-average anomaly of long rains, 2003 - 2018 \n")
plt.xticks([])
plt.yticks([])

ax1 = fig.add_subplot(3, 1, 1)
ax1.set_title("evi average anomaly")
ax1.grid(b=True)
ax1.plot(Years,long_evi_anom,label='evi')
ax1.plot(Years,long_pvi_anom,'r--',label='pvi')
ax1.plot(Years,long_rf_anom,'k--',label='rainfall')
ax1.set_xticks(Years)
ax1.set_xticklabels(Ylabel)
ax1.legend()

ax2 = fig.add_subplot(3, 1, 2)
ax2.set_title("\n\ngsif average anomaly")
ax2.grid(b=True)
ax2.plot(Years,long_gsif_anom,label='gosif')
ax2.plot(Years,long_pvi_anom,'r--',label='pvi')
ax2.plot(Years,long_rf_anom,'k--',label='rainfall')
ax2.set_xticks(Years)
ax2.set_xticklabels(Ylabel)
ax2.legend()


ax3 = fig.add_subplot(3, 1, 3)
ax3.set_title("\n\nnsif average anomaly")
ax3.grid(b=True)
ax3.plot(Years,long_nsif_anom,label='newsif')
ax3.plot(Years,long_pvi_anom,'r--',label='pvi')
ax3.plot(Years,long_rf_anom,'k--',label='rainfall')
ax3.set_xticks(Years)
ax3.set_xticklabels(Ylabel)
ax3.legend()

# fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0, hspace =0.3)#调整子图间距
plt.show()

plt.scatter(short_evi_anom,short_rf_anom)
slope, intercept, r_value, p_value, std_err = stats.linregress(short_evi_anom,short_rf_anom)
print(r_value)
plt.scatter(long_nsif_anom,long_pvi_anom)
slope, intercept, r_value, p_value, std_err = stats.linregress(short_evi_anom,short_rf_anom)
print(r_value)
plt.show()


# DataType ="Crop Mask"
#
# fig = plt.figure(figsize=(15, 10))
#
# plt.title("short rains anomaly correlation coefficient " + " For {}\n".format(DataType), fontsize=16)
# plt.xticks([])
# plt.yticks([])
# for m in range(1, 7):
#
#     ax = fig.add_subplot(2, 3, m)
#     if m == 1:
#         titlename = "ShortRains EVI vs RainFall"
#         Var1List, Var2List = short_evi_anom,short_rf_anom
#         xlabel = "EVI"
#
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#     if m == 2:
#         titlename = "ShortRains GOSIF vs RainFall"
#         Var1List, Var2List = short_gsif_anom,short_rf_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "GOSIF"
#     if m == 3:
#         titlename = "ShortRains NewSIF vs RainFall"
#         Var1List, Var2List = short_nsif_anom,short_rf_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "NEWSIF"
#     if m == 4:
#         titlename = "LongRains EVI vs RainFall"
#         Var1List, Var2List = long_evi_anom,long_rf_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "EVI"
#     if m == 5:
#         titlename = "LongRains GOSIF vs RainFall"
#         Var1List, Var2List = long_gsif_anom,long_rf_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "GOSIF"
#     if m == 6:
#         titlename = "LongRains NSIF vs RainFall"
#         Var1List, Var2List = long_nsif_anom,long_rf_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "NEWSIF"
#
#     scatter = ax.scatter(Var1List, Var2List)
#     titles = titlename
#     ax.set_title(titles, fontsize=10)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel("RainFall")
#     if m <= 3:
#         ax.text(-1.5, 1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#         ax.text(-1.5, 0.7, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#     else:
#
#         ax.text(-1.5, 1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#         ax.text(-1.5, 0.7, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#
# fig.tight_layout()  # 调整整体空白
# # path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\cropatthreshold", DataType + Year + "RF70.jpg")
# # plt.savefig(path2)
# # plt.close()
# plt.show()

#
# DataType ="Crop Mask"
#
# fig = plt.figure(figsize=(15, 10))
#
# plt.title("anomaly correlation coefficient " + " For {}\n".format(DataType), fontsize=16)
# plt.xticks([])
# plt.yticks([])
# for m in range(1, 7):
#
#     ax = fig.add_subplot(2, 3, m)
#     if m == 1:
#         titlename = "ShortRains EVI vs PVI"
#         Var1List, Var2List = short_evi_anom,short_pvi_anom
#         xlabel = "EVI"
#
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#     if m == 2:
#         titlename = "ShortRains GOSIF vs PVI"
#         Var1List, Var2List = short_gsif_anom,short_pvi_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "GOSIF"
#     if m == 3:
#         titlename = "ShortRains NewSIF vs PVI"
#         Var1List, Var2List = short_nsif_anom,short_pvi_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "NEWSIF"
#     if m == 4:
#         titlename = "LongRains EVI vs PVI"
#         Var1List, Var2List = long_evi_anom,long_pvi_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "EVI"
#     if m == 5:
#         titlename = "LongRains GOSIF vs PVI"
#         Var1List, Var2List = long_gsif_anom,long_pvi_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "GOSIF"
#     if m == 6:
#         titlename = "LongRains NSIF vs PVI"
#         Var1List, Var2List = long_nsif_anom,long_pvi_anom
#         slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
#         xlabel = "NEWSIF"
#
#     scatter = ax.scatter(Var1List, Var2List)
#     titles = titlename
#     ax.set_title(titles, fontsize=10)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel("RainFall")
#     if m <= 3:
#         ax.text(-1.5, 1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#         ax.text(-1.5, 0.7, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#     else:
#
#         ax.text(-1.5, 1, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
#         ax.text(-1.5, 0.7, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
#
# fig.tight_layout()  # 调整整体空白
# # path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\cropatthreshold", DataType + Year + "RF70.jpg")
# # plt.savefig(path2)
# # plt.close()
# plt.show()