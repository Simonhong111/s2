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

    short_cmrf_name = ["ShortCMRF" + str(year) for year in range(2006, 2019)]
    long_cmrf_name = ["LongCMRF" + str(year) for year in range(2006, 2019)]
    short_cmpvi_name = ["ShortCMPVI" + str(year) for year in range(2006, 2019)]
    long_cmpvi_name = ["LongCMPVI" + str(year) for year in range(2006, 2019)]

    short_chrf_name = ["ShortCHRF" + str(year) for year in range(2006, 2019)]
    long_chrf_name = ["LongCHRF" + str(year) for year in range(2006, 2019)]
    short_chpvi_name = ["ShortCHPVI" + str(year) for year in range(2006, 2019)]
    long_chpvi_name = ["LongCHPVI" + str(year) for year in range(2006, 2019)]


    short_cmrf =data[short_cmrf_name].to_numpy()
    long_cmrf =data[long_cmrf_name].to_numpy()
    short_cmpvi =data[short_cmpvi_name].to_numpy()*(-1)
    long_cmpvi =data[long_cmpvi_name].to_numpy()*(-1)


    short_chrf =data[short_chrf_name].to_numpy()
    long_chrf =data[long_chrf_name].to_numpy()
    short_chpvi =data[short_chpvi_name].to_numpy()*(-1)
    long_chpvi =data[long_chpvi_name].to_numpy()*(-1)



    return short_cmrf,long_cmrf,short_cmpvi,long_cmpvi,short_chrf,long_chrf,short_chpvi,long_chpvi

path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\CMIP\CMIP_PolyGonAgg_Mask50.csv"

short_cmrf,long_cmrf,short_cmpvi,long_cmpvi,short_chrf,long_chrf,short_chpvi,long_chpvi=monthAnalysis(path)

Years = range(2006,2019)
Ylabel = [str(year)[-2:] for year in range(2006,2019)]


short_cmpvi_anom = (short_cmpvi.mean(axis=0) - short_cmpvi.mean(axis=0).mean())/short_cmpvi.mean(axis=0).std()
long_cmpvi_anom = (long_cmpvi.mean(axis=0) - long_cmpvi.mean(axis=0).mean())/long_cmpvi.mean(axis=0).std()
short_cmrf_anom = (short_cmrf.mean(axis=0) - short_cmrf.mean(axis=0).mean())/short_cmrf.mean(axis=0).std()
long_cmrf_anom = (long_cmrf.mean(axis=0) - long_cmrf.mean(axis=0).mean())/long_cmrf.mean(axis=0).std()

short_chpvi_anom = (short_chpvi.mean(axis=0) - short_chpvi.mean(axis=0).mean())/short_chpvi.mean(axis=0).std()
long_chpvi_anom = (long_chpvi.mean(axis=0) - long_chpvi.mean(axis=0).mean())/long_chpvi.mean(axis=0).std()
short_chrf_anom = (short_chrf.mean(axis=0) - short_chrf.mean(axis=0).mean())/short_chrf.mean(axis=0).std()
long_chrf_anom = (long_chrf.mean(axis=0) - long_chrf.mean(axis=0).mean())/long_chrf.mean(axis=0).std()


fig = plt.figure(figsize=(15, 10))
plt.title("area-average anomaly of short rains, 2006 - 2018 \n")
plt.xticks([])


ax = fig.add_subplot(1, 1, 1)
ax.grid(b=True)
ax.plot(Years,short_cmpvi_anom,'r--',label='cmip5 pvi')
ax.plot(Years,short_cmrf_anom,'k--',label='cmip5 rainfall')
ax.plot(Years,short_chpvi_anom,'r',label='chirps pvi')
ax.plot(Years,short_chrf_anom,'k',label='chirps rainfall')

ax.set_xticks(Years)
ax.set_xticklabels(Ylabel)
ax.legend()

# fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0, hspace =0.3)#调整子图间距
fig = plt.figure(figsize=(15, 10))
plt.title("area-average anomaly of long rains, 2006 - 2018 \n")
plt.xticks([])


ax1 = fig.add_subplot(1, 1, 1)
ax1.grid(b=True)

ax1.plot(Years,long_cmpvi_anom,'r--',label='cmip5 pvi')
ax1.plot(Years,long_cmrf_anom,'k--',label='cmip5 rainfall')
ax1.plot(Years,long_chpvi_anom,'r',label='chirps pvi')
ax1.plot(Years,long_chrf_anom,'k',label='chirps rainfall')
ax1.set_xticks(Years)
ax1.set_xticklabels(Ylabel)
ax1.legend()
# fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0, hspace =0.3)#调整子图间距
plt.show()
# plt.show()