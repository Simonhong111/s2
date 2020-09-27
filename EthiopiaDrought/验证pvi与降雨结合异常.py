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
from scipy import signal
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
    short_pvi =data[short_pvi_name].to_numpy()*(1.0)
    long_pvi =data[long_pvi_name].to_numpy()*(1.0)
    short_gsif =data[short_gsif_name].to_numpy()
    long_gsif =data[long_gsif_name].to_numpy()
    short_nsif =data[short_nsif_name].to_numpy()
    long_nsif =data[long_nsif_name].to_numpy()
    short_Fre =data[short_Fre_name].to_numpy()
    long_Fre = data[long_Fre_name].to_numpy()


    return short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif ,short_Fre ,long_Fre

path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask60_Mean.csv"

short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif ,short_Fre ,long_Fre=monthAnalysis(path)

short_pvi = short_rf/(0.5+short_pvi)
long_pvi = long_rf/(0.5+long_pvi)
nshort_nsif = np.zeros_like(short_nsif,dtype=np.float)
nshort_gsif = np.zeros_like(short_nsif,dtype=np.float)
nlong_nsif = np.zeros_like(short_nsif,dtype=np.float)
nlong_gsif = np.zeros_like(short_nsif,dtype=np.float)

for id in range(nshort_nsif.shape[0]):
    nshort_nsif[id,:] = signal.detrend(short_nsif[id,:])
    nshort_gsif[id, :] = signal.detrend(short_gsif[id, :])
    nlong_nsif[id,:] = signal.detrend(long_nsif[id,:])
    nlong_gsif[id, :] = signal.detrend(long_gsif[id, :])
short_nsif = nshort_nsif
long_nsif = nlong_nsif
short_gsif = nshort_gsif
long_gsif = nlong_gsif


Years = range(2007,2019)
Ylabel = [str(year)[-2:] for year in range(2007,2019)]

short_Fre[short_Fre <=1] = 1
short_Fre[short_Fre >=4] = 4
long_Fre[long_Fre <=1] = 1
long_Fre[long_Fre >=4] = 4

ShortFreSet = set(short_Fre.flatten().tolist())
LongFreSet = set(long_Fre.flatten().tolist())


# 绘制短期的平均值异常

# fig = plt.figure(figsize=(15, 15))
# plt.title("area-average anomaly of short rains, 2007 - 2018 \n\n",fontsize=16)
# plt.xticks([])
# plt.yticks([])
#
# ax1 = fig.add_subplot(4, 1, 1)
# ax1.set_title("average anomaly for frequency <=1 \n")
# ax1.grid(b=True)
#
# mask = np.where(short_Fre == 1)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax1.plot(Years,short_evi_anom,label='evi')
# ax1.plot(Years, short_gsif_anom,'r--', label='gsif')
# ax1.plot(Years, short_nsif_anom, 'k--',label='nsif')
# ax1.plot(Years,short_pvi_anom,label='pvi')
# ax1.plot(Years,short_rf_anom,label='rainfall')
# ax1.set_xticks(Years)
# ax1.set_xticklabels(Ylabel)
# # ax1.legend()
#
# ax2 = fig.add_subplot(4, 1, 2)
# ax2.set_title("average anomaly for frequency 2\n")
# ax2.grid(b=True)
#
# mask = np.where(short_Fre == 2)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax2.plot(Years,short_evi_anom,label='evi')
# ax2.plot(Years, short_gsif_anom,'r--', label='gsif')
# ax2.plot(Years, short_nsif_anom, 'k--',label='nsif')
# ax2.plot(Years,short_pvi_anom,label='pvi')
# ax2.plot(Years,short_rf_anom,label='rainfall')
# ax2.set_xticks(Years)
# ax2.set_xticklabels(Ylabel)
# # ax2.legend()
#
# ax3 = fig.add_subplot(4, 1, 3)
# ax3.set_title("average anomaly for frequency 3\n")
# ax3.grid(b=True)
#
# mask = np.where(short_Fre == 3)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax3.plot(Years,short_evi_anom,label='evi')
# ax3.plot(Years, short_gsif_anom,'r--', label='gsif')
# ax3.plot(Years, short_nsif_anom, 'k--',label='nsif')
# ax3.plot(Years,short_pvi_anom,label='pvi')
# ax3.plot(Years,short_rf_anom,label='rainfall')
# ax3.set_xticks(Years)
# ax3.set_xticklabels(Ylabel)
# # ax3.legend()
#
# ax4 = fig.add_subplot(4, 1, 4)
# ax4.set_title("average anomaly for frequency >=4\n")
# ax4.grid(b=True)
#
# mask = np.where(short_Fre == 4)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax4.plot(Years,short_evi_anom,label='evi')
# ax4.plot(Years, short_gsif_anom,'r--', label='gsif')
# ax4.plot(Years, short_nsif_anom, 'k--',label='nsif')
# ax4.plot(Years,short_pvi_anom,label='pvi')
# ax4.plot(Years,short_rf_anom,label='rainfall')
# ax4.set_xticks(Years)
# ax4.set_xticklabels(Ylabel)
# ax4.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
#
# # fig.tight_layout()#调整整体空白
# plt.subplots_adjust(wspace =0, hspace =0.6)#调整子图间距
# plt.show()
#
#
# #长周期异常
# fig = plt.figure(figsize=(15, 15))
# plt.title("area-average anomaly of long rains, 2007 - 2018 \n\n",fontsize=16)
# plt.xticks([])
# plt.yticks([])
#
# ax1 = fig.add_subplot(4, 1, 1)
# ax1.set_title("average anomaly for frequency <=1 \n")
# ax1.grid(b=True)
#
# mask = np.where(long_Fre == 1)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax1.plot(Years,long_evi_anom,label='evi')
# ax1.plot(Years, long_gsif_anom,'r--', label='gsif')
# ax1.plot(Years, long_nsif_anom, 'k--',label='nsif')
# ax1.plot(Years,long_pvi_anom,label='pvi')
# ax1.plot(Years,long_rf_anom,label='rainfall')
# ax1.set_xticks(Years)
# ax1.set_xticklabels(Ylabel)
# # ax1.legend()
# ax2 = fig.add_subplot(4, 1, 2)
# ax2.set_title("average anomaly for frequency 2 \n")
# ax2.grid(b=True)
#
# mask = np.where(long_Fre == 2)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax2.plot(Years,long_evi_anom,label='evi')
# ax2.plot(Years, long_gsif_anom,'r--', label='gsif')
# ax2.plot(Years, long_nsif_anom, 'k--',label='nsif')
# ax2.plot(Years,long_pvi_anom,label='pvi')
# ax2.plot(Years,long_rf_anom,label='rainfall')
# ax2.set_xticks(Years)
# ax2.set_xticklabels(Ylabel)
#
# ax3 = fig.add_subplot(4, 1, 3)
# ax3.set_title("average anomaly for frequency 3 \n")
# ax3.grid(b=True)
#
# mask = np.where(long_Fre == 3)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax3.plot(Years,long_evi_anom,label='evi')
# ax3.plot(Years, long_gsif_anom,'r--', label='gsif')
# ax3.plot(Years, long_nsif_anom, 'k--',label='nsif')
# ax3.plot(Years,long_pvi_anom,label='pvi')
# ax3.plot(Years,long_rf_anom,label='rainfall')
# ax3.set_xticks(Years)
# ax3.set_xticklabels(Ylabel)
#
# ax4 = fig.add_subplot(4, 1, 4)
# ax4.set_title("average anomaly for frequency >=4 \n")
# ax4.grid(b=True)
#
# mask = np.where(long_Fre == 4)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# ax4.plot(Years,long_evi_anom,label='evi')
# ax4.plot(Years, long_gsif_anom,'r--', label='gsif')
# ax4.plot(Years, long_nsif_anom, 'k--',label='nsif')
# ax4.plot(Years,long_pvi_anom,label='pvi')
# ax4.plot(Years,long_rf_anom,label='rainfall')
# ax4.set_xticks(Years)
# ax4.set_xticklabels(Ylabel)
# ax4.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
#
# # fig.tight_layout()#调整整体空白
# plt.subplots_adjust(wspace =0, hspace =0.6)#调整子图间距
# plt.show()
#
# # 绘制相关性系数

fig = plt.figure(figsize=(15, 10))
Yeas = [str(year)[2:] for year in range(2007,2019)]
plt.title("anomaly correlation coefficient for frequency <=1\n\n", fontsize=16)
plt.xticks([])
plt.yticks([])
mFre = 1
mask = np.where(short_Fre == mFre)
short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()


mask = np.where(long_Fre == mFre)

long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()

for m in range(1, 7):

    ax = fig.add_subplot(2, 3, m)
    if m == 1:
        titlename = "ShortRains EVI vs RainFall\n"
        Var1List, Var2List = short_evi_anom,short_rf_anom
        xlabel = "EVI"

        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
    if m == 2:
        titlename = "ShortRains GOSIF vs RainFall\n"
        Var1List, Var2List = short_gsif_anom,short_rf_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "GOSIF"
    if m == 3:
        titlename = "ShortRains NewSIF vs RainFall\n"
        Var1List, Var2List = short_nsif_anom,short_rf_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "NEWSIF"
    if m == 4:
        titlename = "LongRains EVI vs RainFall\n"
        Var1List, Var2List = long_evi_anom,long_rf_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "EVI"
    if m == 5:
        titlename = "LongRains GOSIF vs RainFall\n"
        Var1List, Var2List = long_gsif_anom,long_rf_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "GOSIF"
    if m == 6:
        titlename = "LongRains NSIF vs RainFall\n"
        Var1List, Var2List = long_nsif_anom,long_rf_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "NEWSIF"

    scatter = ax.scatter(Var1List, Var2List)
    titles = titlename
    ax.set_title(titles, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("RainFall")
    if m <= 3:

        ax.text(0.5, -1.5, "R = {}".format(round(r_value, 3)), fontsize=10, color="r")
        ax.text(0.5, -1.8, "P = {}".format(round(p_value, 3)), fontsize=10, color="r")
        for idx in range(12):
            ax.text(Var1List[idx], Var2List[idx]+0.04, Yeas[idx])
    else:

        ax.text(0.5, -1.5, "R = {}".format(round(r_value, 3)), fontsize=10, color="r")
        ax.text(0.5, -1.8, "P = {}".format(round(p_value, 3)), fontsize=10, color="r")


        for idx in range(12):
            ax.text(Var1List[idx], Var2List[idx]+0.04, Yeas[idx])

plt.subplots_adjust(wspace =0.2, hspace =0.3)#调整子图间距
plt.show()

# 绘制pvi



fig = plt.figure(figsize=(15, 10))
Yeas = [str(year)[2:] for year in range(2007,2019)]
plt.title("Rain/(0.5+PVI) anomaly correlation coefficient for Frequency <=1\n\n", fontsize=16)
plt.xticks([])
plt.yticks([])
mFre = 1
mask = np.where(short_Fre == mFre)
short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()


mask = np.where(long_Fre == mFre)

long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()


for m in range(1, 7):

    ax = fig.add_subplot(2, 3, m)
    if m == 1:
        titlename = "ShortRains EVI vs PVIxRain\n"
        Var1List, Var2List = short_evi_anom,short_pvi_anom
        xlabel = "EVI"

        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
    if m == 2:
        titlename = "ShortRains GOSIF vs PVIxRain\n"
        Var1List, Var2List = short_gsif_anom,short_pvi_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "GOSIF"
    if m == 3:
        titlename = "ShortRains NewSIF vs PVIxRain\n"
        Var1List, Var2List = short_nsif_anom,short_pvi_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "NEWSIF"
    if m == 4:
        titlename = "LongRains EVI vs PVIxRain\n"
        Var1List, Var2List = long_evi_anom,long_pvi_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "EVI"
    if m == 5:
        titlename = "LongRains GOSIF vs PVIxRain\n"
        Var1List, Var2List = long_gsif_anom,long_pvi_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "GOSIF"
    if m == 6:
        titlename = "LongRains NSIF vs PVIxRain\n"
        Var1List, Var2List = long_nsif_anom,long_pvi_anom
        slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
        xlabel = "NEWSIF"

    scatter = ax.scatter(Var1List, Var2List)
    titles = titlename
    ax.set_title(titles, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Rain/(0.5+PVI)")
    if m <= 3:

        ax.text(0.5, -1.1, "R = {}".format(round(r_value, 3)), fontsize=10, color="r")
        ax.text(0.5, -1.3, "P = {}".format(round(p_value, 3)), fontsize=10, color="r")
        for idx in range(12):
            ax.text(Var1List[idx], Var2List[idx]+0.04, Yeas[idx])
    else:

        ax.text(0.5, -1.1, "R = {}".format(round(r_value, 3)), fontsize=10, color="r")
        ax.text(0.5, -1.3, "P = {}".format(round(p_value, 3)), fontsize=10, color="r")


        for idx in range(12):
            ax.text(Var1List[idx], Var2List[idx]+0.04, Yeas[idx])
#
plt.subplots_adjust(wspace =0.2, hspace =0.3)#调整子图间距
plt.show()


# evi_s = []
# evi_long = []
# nsif_s = []
# nsif_long =[]
# gsif_s = []
# gsif_long = []
#
#
# mFre = 2
#
# mask = np.where(short_Fre == mFre)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
#
#
# mask = np.where(long_Fre == mFre)
#
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# Var1List, Var2List = short_evi_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# evi_s.append(r_value)
#
#
# Var1List, Var2List = short_gsif_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# gsif_s.append(r_value)
#
# Var1List, Var2List = short_nsif_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# nsif_s.append(r_value)
#
# Var1List, Var2List = long_evi_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# evi_long.append(r_value)
#
#
# Var1List, Var2List = long_gsif_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# gsif_long.append(r_value)
#
# Var1List, Var2List = long_nsif_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# nsif_s.append(r_value)
#
#
#
# # pvi
# evi_s2 = []
# evi_long2 = []
# nsif_s2 = []
# nsif_long2 =[]
# gsif_s2 = []
# gsif_long2 = []
#
# short_pvi = short_rf/(1+short_pvi)
# long_pvi = long_rf/(1+long_pvi)
# mask = np.where(short_Fre == mFre)
# short_evi_anom = (short_evi[mask[0],:].mean(axis=0) - short_evi[mask[0],:].mean(axis=0).mean()) / short_evi[mask[0],:].mean(axis=0).std()
# short_pvi_anom = (short_pvi[mask[0],:].mean(axis=0) - short_pvi[mask[0],:].mean(axis=0).mean()) / short_pvi[mask[0],:].mean(axis=0).std()
# short_gsif_anom = (short_gsif[mask[0],:].mean(axis=0) - short_gsif[mask[0],:].mean(axis=0).mean()) / short_gsif[mask[0],:].mean(axis=0).std()
# short_nsif_anom = (short_nsif[mask[0],:].mean(axis=0) - short_nsif[mask[0],:].mean(axis=0).mean()) / short_nsif[mask[0],:].mean(axis=0).std()
# short_rf_anom = (short_rf[mask[0],:].mean(axis=0) - short_rf[mask[0],:].mean(axis=0).mean()) / short_rf[mask[0],:].mean(axis=0).std()
#
#
# mask = np.where(long_Fre == mFre)
#
# long_evi_anom = (long_evi[mask[0],:].mean(axis=0) - long_evi[mask[0],:].mean(axis=0).mean()) / long_evi[mask[0],:].mean(axis=0).std()
# long_pvi_anom = (long_pvi[mask[0],:].mean(axis=0) - long_pvi[mask[0],:].mean(axis=0).mean()) / long_pvi[mask[0],:].mean(axis=0).std()
# long_gsif_anom = (long_gsif[mask[0],:].mean(axis=0) - long_gsif[mask[0],:].mean(axis=0).mean()) / long_gsif[mask[0],:].mean(axis=0).std()
# long_nsif_anom = (long_nsif[mask[0],:].mean(axis=0) - long_nsif[mask[0],:].mean(axis=0).mean()) / long_nsif[mask[0],:].mean(axis=0).std()
# long_rf_anom = (long_rf[mask[0],:].mean(axis=0) - long_rf[mask[0],:].mean(axis=0).mean()) / long_rf[mask[0],:].mean(axis=0).std()
#
# Var1List, Var2List = short_evi_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# evi_s2.append(r_value)
#
#
# Var1List, Var2List = short_gsif_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# gsif_s2.append(r_value)
#
# Var1List, Var2List = short_nsif_anom, short_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# nsif_s2.append(r_value)
#
# Var1List, Var2List = long_evi_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# evi_long2.append(r_value)
#
#
# Var1List, Var2List = long_gsif_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# gsif_long2.append(r_value)
#
# Var1List, Var2List = long_nsif_anom, long_rf_anom
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
# nsif_s2.append(r_value)
#
# print(evi_s)
# plt.plot(range(len(evi_s)),evi_s)
# plt.plot(range(len(evi_s2)),evi_s2)
# plt.show()

