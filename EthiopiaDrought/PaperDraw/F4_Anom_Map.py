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
        Evi_path = os.path.join(EVIdir,"{}.{}.tif".format(y,SeasonType))
        assert Evi_path, "The {} does not exist".format(Evi_path)
        Evi_data = gdal.Open(Evi_path).ReadAsArray()

        YearIMGs[:, :, y-st_y] = Evi_data
        YearMask[Evi_data == -9999] = -9999

    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] - Average)/STD
    AnomMap[YearMask == -9999] = -9999

    return AnomMap

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
        SIF_path = os.path.join(SIFdir,"SIF005_{}_{}.nc.tif".format(y,SeasonType))
        # print("0**",SIF_path)
        assert SIF_path, "The {} does not exist".format(SIF_path)
        SIF_data = gdal.Open(SIF_path).ReadAsArray()
        # print(SIF_data[126, 109])
        YearIMGs[:, :, y-st_y] = SIF_data
        YearMask[SIF_data == -9999] = -9999



    # YearIMGsflatten = YearIMGs.reshape(RefW * RefH, YearNum)
    # results = map(signal.detrend, YearIMGsflatten)
    # YearIMGs = np.array(list(results)).reshape(RefH, RefW, YearNum)
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

    RefPath = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D\long_pvi_2007.tif"
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
SeasonType = "Long"
SeasonName = "bega"
Year = 2017
startYear = 2007
P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
              Year=Year,
              SeasonType=SeasonType,startYear=startYear)
RefPath_P = r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage\chirps-v2.0.2007.Long.tif"
RefImg_P = gdal.Open(RefPath_P)
RefGeoTrans_P = RefImg_P.GetGeoTransform()
PX,PY = boundary(RefGeoTrans_P)
P_Anom[P_Anom==-9999] = np.nan
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(a)'.format(SeasonName),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
cmap = plt.get_cmap("RdBu")
cax = plt.imshow(P_Anom,cmap=cmap,vmin=-2,vmax=2)
plt.colorbar(cax,fraction=0.036,pad=0.04)
for i in range(len(PX)):
    plt.plot(PX[i], PY[i], "k")
plt.plot(RX,RY)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_4a{}.png'.format(SeasonName),bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# plt.close(fig)

# Evi
E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
                    Year=Year,
                    SeasonType=SeasonType,startYear=startYear)
RefPath_E =  r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage\2007.Long.tif"
RefImg_E = gdal.Open(RefPath_E)
RefGeoTrans_E = RefImg_E.GetGeoTransform()
EX,EY = boundary(RefGeoTrans_E)
E_Anom[E_Anom==-9999] = np.nan
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(b)'.format(SeasonName),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
Ecmap = plt.get_cmap("RdBu")
Ecax = plt.imshow(E_Anom,cmap=Ecmap,vmin=-2,vmax=2)
plt.colorbar(Ecax,fraction=0.036,pad=0.04)
for i in range(len(EX)):
    plt.plot(EX[i], EY[i], "k")
plt.plot(RX,RY)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_4b{}.png'.format(SeasonName),bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# plt.close(fig)
# SIF
S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage",
                    Year=Year,
                    SeasonType=SeasonType,startYear=startYear)
RefPath_S =  r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\S_SeasonlyAverageImage\SIF005_2007_Long.nc.tif"
RefImg_S = gdal.Open(RefPath_S)
RefGeoTrans_S = RefImg_S.GetGeoTransform()
SX,SY = boundary(RefGeoTrans_S)
S_Anom[S_Anom==-9999] = np.nan
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(c)'.format(SeasonName),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
Scmap = plt.get_cmap("RdBu")
Scax = plt.imshow(S_Anom,cmap=Scmap,vmin=-2,vmax=2)
plt.colorbar(Scax,fraction=0.036,pad=0.04)
for i in range(len(SX)):
    plt.plot(SX[i], SY[i], "k")
plt.plot(RX,RY)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_4c{}.png'.format(SeasonName),bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# plt.close(fig)
# PVI
PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                     Year=Year,
                     SeasonType=SeasonType,startYear=startYear)

RefPath_PV =  r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D\long_pvi_2003.tif"
RefImg_PV = gdal.Open(RefPath_PV)
RefGeoTrans_PV = RefImg_PV.GetGeoTransform()
PVX,PVY = boundary(RefGeoTrans_PV)
PV_Anom[PV_Anom==-9999] = np.nan
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(d)'.format(SeasonName),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
PVcmap = plt.get_cmap("RdBu")
PVcax = plt.imshow(PV_Anom,cmap=PVcmap,vmin=-2,vmax=2)
plt.colorbar(PVcax,fraction=0.036,pad=0.04)
for i in range(len(PVX)):
    plt.plot(PVX[i], PVY[i], "k")
plt.plot(RX,RY)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_4d{}.png'.format(SeasonName),bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
# plt.close(fig)
plt.show()