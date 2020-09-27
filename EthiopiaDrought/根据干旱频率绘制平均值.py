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

    short_evi_name = ["Short"+"EVI"+str(year) for year in range(2003,2019)]
    long_evi_name = ["Long" + "EVI" + str(year) for year in range(2003, 2019)]
    short_rf_name = ["Short" + "RF" + str(year) for year in range(2003, 2019)]
    long_rf_name = ["Long" + "RF" + str(year) for year in range(2003, 2019)]
    short_pvi_name = ["Short" + "PVI" + str(year) for year in range(2003, 2019)]
    long_pvi_name = ["Long" + "PVI" + str(year) for year in range(2003, 2019)]
    short_gsif_name = ["Short" + "GSIF" + str(year) for year in range(2003, 2019)]
    long_gsif_name = ["Long" + "GSIF" + str(year) for year in range(2003, 2019)]
    short_nsif_name = ["Short" + "NSIF" + str(year) for year in range(2003, 2019)]
    long_nsif_name = ["Long" + "NSIF" + str(year) for year in range(2003, 2019)]
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

path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\new\PolyGonAgg_Mask50_Mean.csv"

short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif ,short_Fre ,long_Fre=monthAnalysis(path)
Years = range(2003,2019)
Ylabel = [str(year)[-2:] for year in range(2003,2019)]
ShortFreSet = set(short_Fre.flatten().tolist())
LongFreSet = set(long_Fre.flatten().tolist())
print(ShortFreSet)
print(LongFreSet)

fig = plt.figure(figsize=(15, 10))
plt.title("area-average of short rains, 2003 - 2018 \n")
plt.xticks([])
plt.yticks([])
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_title("evi average")
ax1.grid(b=True)
ax1.plot(Years,short_evi.mean(axis=0))
ax1.set_xticks(Years)
ax1.set_xticklabels(Ylabel)

ax2 = fig.add_subplot(2, 3, 2)
ax2.set_title("pvi average")
ax2.grid(b=True)
ax2.plot(Years,short_pvi.mean(axis=0))
ax2.set_xticks(Years)
ax2.set_xticklabels(Ylabel)

ax3 = fig.add_subplot(2, 3, 3)
ax3.set_title("gsif average")
ax3.grid(b=True)
ax3.plot(Years,short_gsif.mean(axis=0))
ax3.set_xticks(Years)
ax3.set_xticklabels(Ylabel)

ax4 = fig.add_subplot(2, 3, 4)
ax4.set_title("nsif average")
ax4.grid(b=True)
ax4.plot(Years,short_nsif.mean(axis=0))
ax4.set_xticks(Years)
ax4.set_xticklabels(Ylabel)

ax5 = fig.add_subplot(2, 3, 5)
ax5.set_title("rainfall average")
ax5.grid(b=True)
ax5.plot(Years,short_rf.mean(axis=0))
ax5.set_xticks(Years)
ax5.set_xticklabels(Ylabel)
plt.show()
# fig = plt.figure(figsize=(15, 10))
# plt.title("area-average of long rains, 2003 - 2018 \n")
# plt.xticks([])
# plt.yticks([])
# ax1 = fig.add_subplot(2, 3, 1)
# ax1.set_title("evi average")
# ax1.grid(b=True)
# ax1.plot(Years,long_evi.mean(axis=0))
# ax1.set_xticks(Years)
# ax1.set_xticklabels(Ylabel)
#
# ax2 = fig.add_subplot(2, 3, 2)
# ax2.set_title("pvi average")
# ax2.grid(b=True)
# ax2.plot(Years,long_pvi.mean(axis=0))
# ax2.set_xticks(Years)
# ax2.set_xticklabels(Ylabel)
#
# ax3 = fig.add_subplot(2, 3, 3)
# ax3.set_title("gsif average")
# ax3.grid(b=True)
# ax3.plot(Years,long_gsif.mean(axis=0))
# ax3.set_xticks(Years)
# ax3.set_xticklabels(Ylabel)
#
# ax4 = fig.add_subplot(2, 3, 4)
# ax4.set_title("nsif average")
# ax4.grid(b=True)
# ax4.plot(Years,long_nsif.mean(axis=0))
# ax4.set_xticks(Years)
# ax4.set_xticklabels(Ylabel)
#
# ax5 = fig.add_subplot(2, 3, 5)
# ax5.set_title("rainfall average")
# ax5.grid(b=True)
# ax5.plot(Years,long_rf.mean(axis=0))
# ax5.set_xticks(Years)
# ax5.set_xticklabels(Ylabel)
# plt.show()

#
# fig = plt.figure(figsize=(15, 10))
# plt.title("area-average of short rains, 2003 - 2018 \n")
# plt.xticks([])
# plt.yticks([])
# ax1 = fig.add_subplot(2, 3, 1)
# ax1.set_title("evi average")
# ax1.grid(b=True)
#
# for i in ShortFreSet:
#     mask = np.where(short_Fre == i)
#     print(mask[0].shape)
#     ax1.plot(Years,short_evi[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax1.set_xticks(Years)
# ax1.set_xticklabels(Ylabel)
#
# ax2 = fig.add_subplot(2, 3, 2)
# ax2.set_title("pvi average")
# ax2.grid(b=True)
#
# for i in ShortFreSet:
#     mask = np.where(short_Fre == i)
#     print(mask[0].shape)
#     ax2.plot(Years,short_pvi[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax2.set_xticks(Years)
# ax2.set_xticklabels(Ylabel)
#
# ax3 = fig.add_subplot(2, 3, 3)
# ax3.set_title("gsif average")
# ax3.grid(b=True)
#
# for i in ShortFreSet:
#     mask = np.where(short_Fre == i)
#     print(mask[0].shape)
#     ax3.plot(Years,short_gsif[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax3.set_xticks(Years)
# ax3.set_xticklabels(Ylabel)
#
# ax4 = fig.add_subplot(2, 3, 4)
# ax4.set_title("nsif average")
# ax4.grid(b=True)
#
# for i in ShortFreSet:
#     mask = np.where(short_Fre == i)
#     print(mask[0].shape)
#     ax4.plot(Years,short_nsif[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax4.set_xticks(Years)
# ax4.set_xticklabels(Ylabel)
#
# ax5 = fig.add_subplot(2, 3, 5)
# ax5.set_title("rainfall average")
# ax5.grid(b=True)
#
# for i in ShortFreSet:
#     mask = np.where(short_Fre == i)
#     print(mask[0].shape)
#     ax5.plot(Years,short_rf[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " sample num:"+str(mask[0].shape[0]))
# ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
# ax5.set_xticks(Years)
# ax5.set_xticklabels(Ylabel)
# plt.show()

#
# fig = plt.figure(figsize=(15, 10))
# plt.title("area-average of long rains, 2003 - 2018 \n")
# plt.xticks([])
# plt.yticks([])
# ax1 = fig.add_subplot(2, 3, 1)
# ax1.set_title("evi average")
#
# ax1.grid(b=True)
# for i in LongFreSet:
#     mask = np.where(long_Fre == i)
#     print(mask[0].shape)
#     ax1.plot(Years,long_evi[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax1.set_xticks(Years)
# ax1.set_xticklabels(Ylabel)
#
# ax2 = fig.add_subplot(2, 3, 2)
# ax2.set_title("pvi average")
# ax2.grid(b=True)
# for i in LongFreSet:
#     mask = np.where(long_Fre == i)
#     print(mask[0].shape)
#     ax2.plot(Years,long_pvi[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax2.set_xticks(Years)
# ax2.set_xticklabels(Ylabel)
#
# ax3 = fig.add_subplot(2, 3, 3)
# ax3.set_title("gsif average")
# ax3.grid(b=True)
# for i in LongFreSet:
#     mask = np.where(long_Fre == i)
#     print(mask[0].shape)
#     ax3.plot(Years,long_gsif[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax3.set_xticks(Years)
# ax3.set_xticklabels(Ylabel)
#
# ax4 = fig.add_subplot(2, 3, 4)
# ax4.set_title("nsif average")
# ax4.grid(b=True)
# for i in LongFreSet:
#     mask = np.where(long_Fre == i)
#     print(mask[0].shape)
#     ax4.plot(Years,long_nsif[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " num:"+str(mask[0].shape[0]))
#
# ax4.set_xticks(Years)
# ax4.set_xticklabels(Ylabel)
#
# ax5 = fig.add_subplot(2, 3, 5)
# ax5.set_title("rainfall average")
# ax5.grid(b=True)
# for i in LongFreSet:
#     mask = np.where(long_Fre == i)
#     print(mask[0].shape)
#     ax5.plot(Years,long_rf[mask[0],:].mean(axis=0),label="drought fre:"+str(i)+ " sample num:"+str(mask[0].shape[0]))
# ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
# ax5.set_xticks(Years)
# ax5.set_xticklabels(Ylabel)
# plt.show()






