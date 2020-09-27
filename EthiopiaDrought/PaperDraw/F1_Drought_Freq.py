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

def AnomalyFrequency(chirpsdirectory,MonthType):
    start = datetime.strptime("1981-01-01", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    refpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage\chirps-v2.0.1981.Long.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    Width, Height = reference.RasterXSize, reference.RasterYSize
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}.{}.tif".format(str(dt.year),MonthType))

        if os.path.exists(chirps_file):
            bandNum += 1
        assert os.path.exists(chirps_file),"this file does not exist：".format(chirps_file)

    frequency = np.zeros((Height, Width))
    maskarr = np.zeros((Height, Width))
    multidarr = np.zeros((Height, Width, bandNum))

    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}.{}.tif".format(str(dt.year),MonthType))

        print(dt, chirps_file)
        chirps = gdal.Open(chirps_file).ReadAsArray()

        multidarr[:, :, band_id] = chirps
        maskarr[chirps == -9999] = -9999
        band_id += 1

    stdMatrix = multidarr.std(axis=2)
    averageMatrix = multidarr.mean(axis=2)

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"chirps-v2.0.{}.{}.tif".format(str(dt.year),MonthType))

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
    multidarr[maskarr == -9999] = -9999

    return frequency,multidarr


ShortFrequency,ShortPmat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage","Short")
LongFrequency,LongPmat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage","Long")

# 绘图
NW = [68,42,132,157]

vmin = ShortFrequency[ShortFrequency > 0].min()
vmax = ShortFrequency[ShortFrequency > 0].max()
print("maxmin value", vmin, vmax)

vmin = LongFrequency[LongFrequency > 0].min()
vmax = LongFrequency[LongFrequency > 0].max()
print("long maxmin value", vmin, vmax)

ShortFrequency[ShortFrequency ==-9999] = np.nan
LongFrequency[LongFrequency==-9999] = np.nan


# Short Rains
fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(1,2,1)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
cmap = plt.get_cmap("coolwarm")
boundaries = [1.5,2.5, 3.5, 4.5, 5.5,6.5,7.5,8.5,9.5,10.5]
norm = mpl.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
cax = ax.imshow(ShortFrequency, cmap=cmap, norm=norm)
cmap.set_under('white')
cbar = plt.colorbar(cax, ax=ax,fraction=0.0355, pad=0.04)
cbar.ax.get_yaxis().set_ticks([])
cbar.ax.get_yaxis().set_ticks([2,3,4,5,6,7,8,9,10])
cbar.ax.get_yaxis().set_ticklabels([2,3,4,5,6,7,8,9,10])
cbar.ax.get_yaxis().set_tick_params(direction='in')
ax.text(0.1, 0.95,'a) Bega',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)

# Long Rains
ax2 = fig.add_subplot(1,2,2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axis('off')
cax2 = ax2.imshow(LongFrequency, cmap=cmap, norm=norm)
cmap.set_under('white')
cbar2 = plt.colorbar(cax2,ax=ax2,fraction=0.0355, pad=0.04)
cbar2.ax.get_yaxis().set_ticks([])
cbar2.ax.get_yaxis().set_ticks([2,3,4,5,6,7,8,9,10])
cbar2.ax.get_yaxis().set_ticklabels([2,3,4,5,6,7,8,9,10])
cbar2.ax.get_yaxis().set_tick_params(direction='in')
plt.text(0.1, 0.95,'b) Kremt',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.05)

plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_1_Long.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.0)
plt.show()







