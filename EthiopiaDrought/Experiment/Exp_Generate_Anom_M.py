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
from Experiment.Exp_Write_Image import write_Img as Write
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

def SIFAnomMap(SIFdir,OutDir,Year,Month,startYear=2003):
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
        SIF_path = os.path.join(SIFdir,"SIF005_{}{}.nc.tif".format(y,Month))
        # print("0**",SIF_path)
        assert SIF_path, "The {} does not exist".format(SIF_path)
        SIF_data = gdal.Open(SIF_path).ReadAsArray()
        # print(SIF_data[126, 109])
        YearIMGs[:, :, y-st_y] = SIF_data
        YearMask[SIF_data == -9999] = -9999

    Average = YearIMGs.mean(axis=2)
    STD = YearIMGs.std(axis=2)
    AnomMap = (YearIMGs[:,:,Year-st_y] -Average)/STD
    AnomMap[YearMask == -9999] = -9999
    print(AnomMap[126, 109])
    out_path = os.path.join(OutDir, "SIF005_{}{}.nc.tif".format(Year, Month))
    Write(AnomMap, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)


for year in range(2007,2019):
    for m in ["02","03","04","05","06", "07", "08", "09"]:
        SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_SIF_V3_Ethiopia",
                   r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\MDetrend_SIF_V3_Ethiopia",
                   year, m, startYear=2007)





