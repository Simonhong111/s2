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
def CHIPRSAnomMap(Chirpsdir,Year,COL,ROW,SeasonType='05',startYear=2003):
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
        YearIMGs[:, :, y-st_y] = ch_data
        YearMask[ch_data == -9999] = -9999


    valuelist = YearIMGs[ROW,COL,:]
    print("P"+SeasonType,valuelist)
    valuelist = (valuelist - valuelist.mean())/valuelist.std()
    value = valuelist[Year-startYear]


    return value,valuelist

def EVIAnomMap(EVIdir,Year,COL,ROW,SeasonType='Short',startYear=2003):
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

        YearIMGs[:, :, y - st_y] = Evi_data
        YearMask[Evi_data == -9999] = -9999

    valuelist = YearIMGs[ROW, COL, :]
    valuelist = (valuelist - valuelist.mean()) / valuelist.std()
    value = valuelist[Year - startYear]

    return value,valuelist

def SIFAnomMap(SIFdir,Year,COL,ROW,SeasonType='Short',startYear=2003):
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
        YearIMGs[:, :, y - st_y] = SIF_data
        YearMask[SIF_data == -9999] = -9999

    valuelist = YearIMGs[ROW, COL, :]
    valuelist = np.array(list(signal.detrend(valuelist)))
    valuelist = (valuelist - valuelist.mean()) / valuelist.std()
    value = valuelist[Year - startYear]





    return value,valuelist


def PVIAnomMap(PVIdir,Year,COL,ROW,SeasonType='Short',startYear=2003):
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
        # print(PVI_data[70,89])
        YearIMGs[:, :, y - st_y] = PVI_data
        YearMask[PVI_data == -9999] = -9999

    valuelist = YearIMGs[ROW, COL, :]
    valuelist = (valuelist - valuelist.mean()) / valuelist.std()
    value = valuelist[Year - startYear]

    # breakpoint()
    return value,valuelist


COL = 76
ROW = 103
Year = 2018
startYear = 2007
Months= ['06','07','08','09']
# Months= ['02','03','04','05']
FM = []
FML = []
SM =[]
SML = []
TML =[]
TM = []
FRML = []
FRM =[]
for m in Months:
    ch_v,ch_list = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\Chirps2", Year, COL, ROW, SeasonType=m, startYear=startYear)
    e_v,e_list = EVIAnomMap(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia", Year, COL, ROW, SeasonType=m, startYear=startYear)
    s_v,s_list = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia", Year, COL, ROW,
               SeasonType=m, startYear=startYear)
    p_v,p_list = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D", Year, COL, ROW, SeasonType='Long',
               startYear=startYear)
    FM.append(ch_v)
    FML.append(ch_list)
    SM.append(e_v)
    SML.append(e_list)
    TM.append(s_v)
    TML.append(s_list)
    FRM.append(p_v)
    FRML.append(p_list)
fig = plt.figure(1)
plt.title("1")
plt.grid(True)
plt.scatter([Year],[FM[0]])
plt.plot([[i+startYear] for i in range(2019-startYear)],FML[0],label="1")
plt.plot([[i+startYear] for i in range(2019-startYear)],FML[1],label="2")
plt.plot([[i+startYear] for i in range(2019-startYear)],FML[2],label='3')
plt.plot([[i+startYear] for i in range(2019-startYear)],FML[3],label='4')


plt.legend()
fig = plt.figure(2)
plt.title("2")
plt.grid(True)
plt.scatter([Year],[SM[1]])
plt.plot([[i+startYear] for i in range(2019-startYear)],SML[0],label="1")
plt.plot([[i+startYear] for i in range(2019-startYear)],SML[1],label="2")
plt.plot([[i+startYear] for i in range(2019-startYear)],SML[2],label="3")
plt.plot([[i+startYear] for i in range(2019-startYear)],SML[3],label="4")
plt.legend()

fig = plt.figure(3)
plt.title("3")
plt.grid(True)
plt.scatter([Year],[TM[2]])
plt.plot([[i+startYear] for i in range(2019-startYear)],TML[0],label="1")
plt.plot([[i+startYear] for i in range(2019-startYear)],TML[1],label="2")
plt.plot([[i+startYear] for i in range(2019-startYear)],TML[2],label="3")
plt.plot([[i+startYear] for i in range(2019-startYear)],TML[3],label="4")
# TML = np.array(TML).mean(axis=0)
# print(TML.shape)
# TML = (TML - TML.mean())/TML.std()
# plt.plot([[i+startYear] for i in range(2019-startYear)],TML,label="5")

plt.legend()
fig = plt.figure(4)
plt.title("4")
plt.grid(True)
plt.scatter([Year],[FRM[3]])
plt.plot([[i+startYear] for i in range(2019-startYear)],FRML[0],label="1")
plt.plot([[i+startYear] for i in range(2019-startYear)],FRML[1],label="2")
plt.plot([[i+startYear] for i in range(2019-startYear)],FRML[2],label="3")
plt.plot([[i+startYear] for i in range(2019-startYear)],FRML[3],label="4")
plt.legend()
plt.show()



