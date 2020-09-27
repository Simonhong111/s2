from osgeo import gdal,osr,ogr
import os,glob
import numpy as np
from dateutil import rrule
from datetime import *
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import signal
from scipy import signal
from Experiment.Exp_Ethiopia_SubRegion_Shape import boundary
import matplotlib.gridspec as gridspec
import ctypes
ctypes.WinDLL("kernel32.dll")
ctypes.WinDLL("msvcrt.dll")
ctypes.WinDLL("user32.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libwinpthread-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgcc_s_seh-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgomp-1.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\cblas\.libs\libgslcblas-0.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\.libs\libgsl-25.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\libeemd\libeemd.so")
from pyeemd import ceemdan,eemd
from pyeemd.utils import plot_imfs
import matplotlib.pyplot as plt
def CHIPRSAnomMap(Chirpsdir,Year,SeasonType='Short',startYear=2003):
    """
    :param Chirpsdir: Seasonal average image dir
    :param Year: the year of your anomaly map
    :param SeasonType: Short Rains or Long Rains
    :return: Anomaly Map in the Year you given
    """

    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage\chirps-v2.0.2007.Long.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW, RefH = RefImg.RasterXSize, RefImg.RasterYSize

    st_y = startYear
    en_y = 2018
    YearNum = en_y - st_y + 1
    YearIMGs = np.zeros(shape=(RefH, RefW, YearNum), dtype=np.float)
    YearMask = np.zeros(shape=(RefH, RefW))

    # data cover 2003-2018
    for y in range(st_y, en_y+1):
        ch_path = os.path.join(Chirpsdir,"chirps-v2.0.{}.{}.tif".format(y,SeasonType))
        assert ch_path, "The {} does not exist".format(ch_path)

        ch_data = gdal.Open(ch_path).ReadAsArray()
        # print(ch_data[126,109])
        YearIMGs[:, :, y-st_y] = ch_data
        YearMask[ch_data == -9999] = -9999
    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] - Average)/STD
    AnomMap[YearMask == -9999] = -9999
    # print(AnomMap[126,109])
    return AnomMap

def EVIAnomMap(EVIdir,Year,SeasonType='Short',startYear=2003):
    """
    :param Chirpsdir: Seasonal average image dir
    :param Year: the year of your anomaly map
    :param SeasonType: Short Rains or Long Rains
    :return: Anomaly Map in the Year you given
    """

    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage\2007.Long.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW, RefH = RefImg.RasterXSize, RefImg.RasterYSize
    st_y = startYear
    en_y = 2018
    YearNum = en_y - st_y + 1

    YearIMGs = np.zeros(shape=(RefH, RefW, YearNum), dtype=np.float)
    YearMask = np.zeros(shape=(RefH, RefW))

    # data cover 2003-2018
    for y in range(st_y, en_y+1):
        Evi_path = os.path.join(EVIdir,"{}.{}.01.tif".format(y,SeasonType))
        assert Evi_path, "The {} does not exist".format(Evi_path)
        Evi_data = gdal.Open(Evi_path).ReadAsArray()

        YearIMGs[:, :, y-st_y] = Evi_data
        YearMask[Evi_data == -9999] = -9999

    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] - Average)/STD
    AnomMap[YearMask == -9999] = -9999

    return AnomMap
def fit(y):
    imfs = ceemdan(y,S_number=4,num_siftings=50)
    trend = imfs[-1]+imfs[-2]
    dtrend = y-trend
    return dtrend
def SIFAnomMap(SIFdir,Year,SeasonType='Short',startYear=2003):
    """
    :param Chirpsdir: Seasonal average image dir
    :param Year: the year of your anomaly map
    :param SeasonType: Short Rains or Long Rains
    :return: Anomaly Map in the Year you given
    """

    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage\SIF005_2007_Long.nc.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW, RefH = RefImg.RasterXSize, RefImg.RasterYSize
    st_y = startYear
    en_y = 2018
    YearNum = en_y - st_y + 1

    YearIMGs = np.zeros(shape=(RefH, RefW, YearNum), dtype=np.float)
    YearMask = np.zeros(shape=(RefH, RefW))

    # data cover 2003-2018
    for y in range(st_y, en_y+1):
        SIF_path = os.path.join(SIFdir,"SIF005_{}{}.nc.tif".format(y,SeasonType))
        # print("0**",SIF_path)
        assert SIF_path, "The {} does not exist".format(SIF_path)
        SIF_data = gdal.Open(SIF_path).ReadAsArray()
        # print(SIF_data[126, 109])
        YearIMGs[:, :, y-st_y] = SIF_data
        YearMask[SIF_data == -9999] = -9999



    YearIMGsflatten = YearIMGs.reshape(RefW * RefH, YearNum)
    results = map(signal.detrend, YearIMGsflatten)
    YearIMGs = np.array(list(results)).reshape(RefH, RefW, YearNum)
    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] -Average)/STD
    AnomMap[YearMask == -9999] = -9999
    print(AnomMap[126, 109])
    return AnomMap


def PVIAnomMap(PVIdir,Year,SeasonType='Short',startYear=2003):
    """
    :param Chirpsdir: Seasonal average image dir
    :param Year: the year of your anomaly map
    :param SeasonType: Short Rains or Long Rains
    :return: Anomaly Map in the Year you given
    """

    RefPath = r"D:\Cornell\EthiopianDrought\AData\PVI10day\long_pvi_2003.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW, RefH = RefImg.RasterXSize, RefImg.RasterYSize
    st_y = startYear
    en_y = 2018
    YearNum = en_y - st_y + 1

    YearIMGs = np.zeros(shape=(RefH, RefW, YearNum), dtype=np.float)
    YearMask = np.zeros(shape=(RefH, RefW))

    # data cover 2003-2018
    for y in range(st_y, en_y+1):
        PVI_path = os.path.join(PVIdir,"{}_pvi_{}.tif".format(SeasonType,y))
        assert PVI_path, "The {} does not exist".format(PVI_path)
        PVI_data = gdal.Open(PVI_path).ReadAsArray()
        # print(PVI_data[126, 109])
        YearIMGs[:, :, y-st_y] = PVI_data*(-1.0)
        YearMask[PVI_data == -9999] = -9999


    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] - Average)/STD
    AnomMap[YearMask == -9999] = -9999

    # print(AnomMap[126, 109])
    return AnomMap


RX = [68,132,132,68,68]
RY = [42,42,157,157,42]



# precipitation
SeasonType = "05"
SeasonName = "Belg"
Year = 2018
startYear = 2007
if SeasonType == "Short":
    titles = ["(a)","(b)","(c)","(d)"]
else:
    titles = ["(e)","(f)","(g)","(h)"]


fig = plt.figure(figsize=(10, 6))
plt.title(str(Year)+SeasonType)
P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\Chirps2",
                       Year=Year,
                       SeasonType=SeasonType, startYear=startYear)
RefPath_P = r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage\chirps-v2.0.2007.Long.tif"
RefImg_P = gdal.Open(RefPath_P)
RefGeoTrans_P = RefImg_P.GetGeoTransform()
PX, PY = boundary(RefGeoTrans_P)

P_Anom[P_Anom == -9999] = np.nan
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.text(0.05, 0.95, '{}'.format(titles[0]), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='left',
         verticalalignment='center',
         transform=ax1.transAxes)
ax1.text(0.95, 0.95, '{} P'.format(SeasonName), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax1.transAxes)
cmap = plt.get_cmap("RdBu")
cax = ax1.imshow(P_Anom, cmap=cmap, vmin=-2, vmax=2)
plt.colorbar(cax, ax=ax1, fraction=0.0362, pad=0.04)
for i in range(len(PX)):
    ax1.plot(PX[i], PY[i], "k")
ax1.plot(RX, RY)

# Evi
E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia",
                    Year=Year,
                    SeasonType=SeasonType, startYear=startYear)
RefPath_E = r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage\2007.Long.tif"
RefImg_E = gdal.Open(RefPath_E)
RefGeoTrans_E = RefImg_E.GetGeoTransform()
EX, EY = boundary(RefGeoTrans_E)
E_Anom[E_Anom == -9999] = np.nan
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.text(0.05, 0.95, '{}'.format(titles[1]), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='left',
         verticalalignment='center',
         transform=ax2.transAxes)
ax2.text(0.95, 0.95, '{} EVI'.format(SeasonName), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax2.transAxes)
Ecmap = plt.get_cmap("RdBu")
Ecax = ax2.imshow(E_Anom, cmap=Ecmap, vmin=-2, vmax=2)
plt.colorbar(Ecax, ax=ax2, fraction=0.0362, pad=0.04)
for i in range(len(EX)):
    ax2.plot(EX[i], EY[i], "k")
ax2.plot(RX, RY)

# SIF
S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia",
                    Year=Year,
                    SeasonType=SeasonType, startYear=startYear)
RefPath_S = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage\SIF005_2007_Long.nc.tif"
RefImg_S = gdal.Open(RefPath_S)
RefGeoTrans_S = RefImg_S.GetGeoTransform()
SX, SY = boundary(RefGeoTrans_S)
S_Anom[S_Anom == -9999] = np.nan
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xticks([])
ax3.set_yticks([])

ax3.text(0.05, 0.95, '{}'.format(titles[2]), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='left',
         verticalalignment='center',
         transform=ax3.transAxes)
ax3.text(0.95, 0.95, '{} SIF'.format(SeasonName), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax3.transAxes)
Scmap = plt.get_cmap("RdBu")
Scax = ax3.imshow(S_Anom, cmap=Scmap, vmin=-2, vmax=2)
plt.colorbar(Scax, ax=ax3, fraction=0.0362, pad=0.04)
for i in range(len(SX)):
    ax3.plot(SX[i], SY[i], "k")
ax3.plot(RX, RY)

# PVI
if SeasonName == "Belg":
    ST = "Short"
if SeasonName =="Kiremit":
    ST = "Long"
PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                     Year=Year,
                     SeasonType=ST, startYear=startYear)

RefPath_PV = r"D:\Cornell\EthiopianDrought\AData\PVI5day\long_pvi_2003.tif"
RefImg_PV = gdal.Open(RefPath_PV)
RefGeoTrans_PV = RefImg_PV.GetGeoTransform()
PVX, PVY = boundary(RefGeoTrans_PV)
PV_Anom[PV_Anom == -9999] = np.nan
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.text(0.05, 0.95, '{}'.format(titles[3]), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='left',
         verticalalignment='center',
         transform=ax4.transAxes)
ax4.text(0.95, 0.95, '{} PVI'.format(SeasonName), fontdict={'fontname': 'Times New Roman', 'fontsize': 16},
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax4.transAxes)
PVcmap = plt.get_cmap("RdBu")
PVcax = ax4.imshow(PV_Anom, cmap=PVcmap, vmin=-2, vmax=2)
plt.colorbar(PVcax, ax=ax4, fraction=0.0362, pad=0.04)
for i in range(len(PVX)):
    ax4.plot(PVX[i], PVY[i], "k")
ax4.plot(RX, RY)


plt.subplots_adjust(left=0.1, right=0.9, bottom=0, top=0.9, wspace=0.15, hspace=0.05)


plt.show()