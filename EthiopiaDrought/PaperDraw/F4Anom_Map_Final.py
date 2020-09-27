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
SeasonType = "Long"
SeasonName = "Kiremit"
# Year = 2011
startYear = 2007
if SeasonType == "Short":
    titles = ["(a)","(b)","(c)","(d)"]
else:
    titles = ["(e)","(f)","(g)","(h)"]

# fig = plt.figure(constrained_layout=True)
fig = plt.figure()
gs = fig.add_gridspec(6, 8)
cmap = plt.get_cmap("RdBu")
for Year in range(5):

        SeasonType = "Short"
        ax1 = fig.add_subplot(gs[Year,0])
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_ylabel("{}".format(Year+startYear), fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
        P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
                               Year=Year+startYear,
                               SeasonType=SeasonType, startYear=startYear)
        P_Anom[P_Anom == -9999] = np.nan
        Pcax = ax1.imshow(P_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # Evi
        ax2 = fig.add_subplot(gs[Year,1])
        ax2.set_xticks([])
        ax2.set_yticks([])
        E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
                            Year=Year+startYear,
                            SeasonType=SeasonType, startYear=startYear)
        E_Anom[E_Anom == -9999] = np.nan
        Ecax = ax2.imshow(E_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # SIF
        ax3 = fig.add_subplot(gs[Year,2])
        ax3.set_xticks([])
        ax3.set_yticks([])
        S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage",
                            Year=Year+startYear,
                            SeasonType=SeasonType, startYear=startYear)
        S_Anom[S_Anom == -9999] = np.nan
        Scax = ax3.imshow(S_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # PVI
        ax4 = fig.add_subplot(gs[Year,3])
        ax4.set_xticks([])
        ax4.set_yticks([])
        PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                             Year=Year+startYear,
                             SeasonType=SeasonType, startYear=startYear)
        PV_Anom[PV_Anom == -9999] = np.nan
        PVIcax = ax4.imshow(PV_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)


        #  *****************
        SeasonType = "Long"
        ax5 = fig.add_subplot(gs[Year,4])
        ax5.set_xticks([])
        ax5.set_yticks([])
        P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
                               Year=Year+startYear,
                               SeasonType=SeasonType, startYear=startYear)
        P_Anom[P_Anom == -9999] = np.nan
        Pcax = ax5.imshow(P_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # Evi
        ax6 = fig.add_subplot(gs[Year,5])
        ax6.set_xticks([])
        ax6.set_yticks([])
        E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
                            Year=Year+startYear,
                            SeasonType=SeasonType, startYear=startYear)
        E_Anom[E_Anom == -9999] = np.nan
        Ecax = ax6.imshow(E_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # SIF
        ax7 = fig.add_subplot(gs[Year,6], sharex=ax1)
        ax7.set_xticks([])
        ax7.set_yticks([])
        S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage",
                            Year=Year+startYear,
                            SeasonType=SeasonType, startYear=startYear)
        S_Anom[S_Anom == -9999] = np.nan
        Scax = ax7.imshow(S_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

        # PVI
        ax8 = fig.add_subplot(gs[Year,7])
        ax8.set_xticks([])
        ax8.set_yticks([])
        PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                             Year=Year+startYear,
                             SeasonType=SeasonType, startYear=startYear)
        PV_Anom[PV_Anom == -9999] = np.nan
        PVIcax = ax8.imshow(PV_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)
# Final Year
SeasonType = "Short"
Year = 5
ax1 = fig.add_subplot(gs[Year,0])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_ylabel("{}".format(Year+startYear), fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
ax1.set_xlabel("P", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
                       Year=Year+startYear,
                       SeasonType=SeasonType, startYear=startYear)
P_Anom[P_Anom == -9999] = np.nan
Pcax = ax1.imshow(P_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# Evi
ax2 = fig.add_subplot(gs[Year,1])
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel("EVI", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
                    Year=Year+startYear,
                    SeasonType=SeasonType, startYear=startYear)
E_Anom[E_Anom == -9999] = np.nan
Ecax = ax2.imshow(E_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# SIF
ax3 = fig.add_subplot(gs[Year,2])
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel("SIF", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage",
                    Year=Year+startYear,
                    SeasonType=SeasonType, startYear=startYear)
S_Anom[S_Anom == -9999] = np.nan
Scax = ax3.imshow(S_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# PVI
ax4 = fig.add_subplot(gs[Year,3])
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlabel("PVI", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                     Year=Year+startYear,
                     SeasonType=SeasonType, startYear=startYear)
PV_Anom[PV_Anom == -9999] = np.nan
PVIcax = ax4.imshow(PV_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)


#  *****************
SeasonType = "Long"
ax5 = fig.add_subplot(gs[Year,4])
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlabel("P", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})

P_Anom = CHIPRSAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage",
                       Year=Year+startYear,
                       SeasonType=SeasonType, startYear=startYear)
P_Anom[P_Anom == -9999] = np.nan
Pcax = ax5.imshow(P_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# Evi
ax6 = fig.add_subplot(gs[Year,5])
ax6.set_xticks([])
ax6.set_yticks([])
ax6.set_xlabel("EVI", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
E_Anom = EVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\Evi_Data\E_SeasonlyAverageImage",
                    Year=Year+startYear,
                    SeasonType=SeasonType, startYear=startYear)
E_Anom[E_Anom == -9999] = np.nan
Ecax = ax6.imshow(E_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# SIF
ax7 = fig.add_subplot(gs[Year,6], sharex=ax1)
ax7.set_xticks([])
ax7.set_yticks([])
ax7.set_xlabel("SIF", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
S_Anom = SIFAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_S_SeasonlyAverageImage",
                    Year=Year+startYear,
                    SeasonType=SeasonType, startYear=startYear)
S_Anom[S_Anom == -9999] = np.nan
Scax = ax7.imshow(S_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

# PVI
ax8 = fig.add_subplot(gs[Year,7])
ax8.set_xticks([])
ax8.set_yticks([])
ax8.set_xlabel("PVI", fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
PV_Anom = PVIAnomMap(r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D",
                     Year=Year+startYear,
                     SeasonType=SeasonType, startYear=startYear)
PV_Anom[PV_Anom == -9999] = np.nan
PVIcax = ax8.imshow(PV_Anom[42:157,68:132], cmap=cmap, vmin=-2, vmax=2)

fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.9, 0.11, 0.01, 0.77])
fig.colorbar(PVIcax, cax=cbar_ax)




# plt.tight_layout()
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0, top=0.9, wspace=0, hspace=0)
fig.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_4.png',
            bbox_inches='tight', dpi=fig.dpi, pad_inches=0.05)
plt.show()





