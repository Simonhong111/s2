from osgeo import gdal,osr,ogr
import os
import numpy as np
from dateutil import rrule
from datetime import *
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as patches
def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

def AnomalyFrequency(chirpsdirectory,Monthtype):
    start = datetime.strptime("1981-01-01", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    refpath = r"D:\Cornell\EthiopianDrought\CHIRPSRainSeasonCom\SeasonCom\chirps-v2.0.Rains_1981.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    Width, Height = reference.RasterXSize, reference.RasterYSize
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}_Rains_{}.tif".format(Monthtype,str(dt.year)))

        if os.path.exists(chirps_file):
            bandNum += 1
        assert os.path.exists(chirps_file),"this file does not exist：".format(chirps_file)

    frequency = np.zeros((Height, Width))
    maskarr = np.zeros((Height, Width))
    multidarr = np.zeros((Height, Width, bandNum))

    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}_Rains_{}.tif".format(Monthtype,str(dt.year)))

        print(dt, chirps_file)
        chirps = gdal.Open(chirps_file).ReadAsArray()

        multidarr[:, :, band_id] = chirps
        maskarr[chirps == -9999] = -9999
        band_id += 1

    stdMatrix = multidarr.std(axis=2)
    averageMatrix = multidarr.mean(axis=2)

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}_Rains_{}.tif".format(Monthtype,str(dt.year)))

        chirps = gdal.Open(chirps_file).ReadAsArray()

        mask2 = (stdMatrix > 0) & (maskarr > -9999) & (chirps > -9999)
        mask2 = np.where(mask2)
        anomalyMap = np.zeros(shape=chirps.shape,dtype=np.float)
        anomalyMap[mask2] = (chirps[mask2] - averageMatrix[mask2]) / stdMatrix[mask2]
        frequency[anomalyMap < -1] +=1

        del chirps
        del mask2
        del anomalyMap
    mask3 = (stdMatrix <= 0) | (maskarr == -9999)
    frequency[mask3] = -9999
    multidarr[maskarr ==-9999] = -9999

    return frequency,multidarr


ShortFrequency,ShortPMat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\CHIRPSRainSeasonCom\SeasonCom","Short")
LongFrequency,LongPMat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\CHIRPSRainSeasonCom\SeasonCom","Long")

# 绘图
NW = [68,42,132,157]

ShortNWsubmat = ShortPMat[NW[1]:NW[3],NW[0]:NW[2],:]
LongNWsubmat = LongPMat[NW[1]:NW[3],NW[0]:NW[2],:]


SNWTS = []
LNWTS = []

H,W ,Years = ShortNWsubmat.shape
for y in range(Years):
    SNWtp = ShortNWsubmat[:,:,y]
    SNWTS.append(SNWtp[SNWtp >-9999].mean())
    LNWtp = LongNWsubmat[:, :, y]
    LNWTS.append(LNWtp[LNWtp > -9999].mean())


fig = plt.figure(figsize=(5, 3))
# plt.xticks([])
# plt.yticks([])
# plt.axis('off')
plt.plot(range(Years),SNWTS,'k-*',label="precipitation")
plt.xticks(range(0,39,5),[str(1981+y)[2:] for y in range(0,39,5)])
plt.text(0.1, 0.95,'c) bega',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
# plt.legend(loc=2)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_1c.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()
fig = plt.figure(figsize=(5, 3))
# plt.xticks([])
# plt.yticks([])
# plt.axis('off')
plt.plot(range(Years),LNWTS,'k-*',label="precipitation")

plt.xticks(range(0,39,5),[str(1981+y)[2:] for y in range(0,39,5)])
# plt.legend(loc=2)
plt.text(0.1, 0.95,'d) kremt',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)

plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_1d.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()



