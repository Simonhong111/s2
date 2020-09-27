from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import signal

dirname = r'D:\Cornell\GOSIFV002Clip'

month = 2

i = 0
Value = []
Tm = []


start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()

for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
    i += 1
    Tm.append(dt.year)
    file = os.path.join(dirname,"GOSIF_{}.M{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
    print(file)
    if dt.year !=2017:
        continue
    sif = gdal.Open(file).ReadAsArray()
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', 'MCD12C1.A{}001.*.tif'.format(str(dt.year))))[0]
    crop = gdal.Open(cropland).ReadAsArray()
    sif[crop != 12] = -9999
    sifarr = sif[(sif > -999) & (sif < 32766)]*0.0001
    plt.hist(sifarr,bins=100)
    # mean = sif[(sif > -999) & (sif < 32766)].mean()
    #
    # print(dt, sif[(sif > -999) & (sif < 32766)].std())
    # Value.append(mean)
#     fig = plt.figure(i)
#     plt.title(str(dt.year)+str(dt.month).zfill(2))
#     cmp = plt.cm.get_cmap('Greys')
#     cmp.set_under('white')
#     plt.imshow(sif,cmap=cmp,vmin=0,vmax=0.8)
#     plt.colorbar()
# plt.show()
def fit(y):
    return signal.detrend(y)

# mean = np.array(Value).mean()
# # Value1 = list(fit(Value[0:6]))
# # Value2 = list(fit(Value[6:]))
# #
# #
# # Value1.extend(Value2)
# # Value1 = np.array(Value1)
# # Value1 = (Value1 - np.mean(Value1))/np.std(Value1)
# Value1 = fit(Value)+mean
# print(Tm)
#
# plt.plot(Tm,Value)
# plt.plot(Tm,Value1,'r')
# plt.scatter(Tm,Value)
# plt.scatter(Tm,Value1)
# # plt.plot(Tm,[0]*len(Tm))
# print(Value)
plt.show()

os.path.ex

