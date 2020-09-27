from osgeo import gdal,osr,ogr
import os
import numpy as np
from dateutil import rrule
from datetime import *
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
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

print("duan",np.where(ShortFrequency >=6)[0].shape)
print("duan",np.where(ShortFrequency != -9999)[0].shape)
print(ShortFrequency.shape[0]*ShortFrequency.shape[1])
print("duan",np.where(LongFrequency >=6)[0].shape)
print("duan",np.where(LongFrequency != -9999)[0].shape)
print(LongFrequency.shape[0]*LongFrequency.shape[1])

# 绘图
NW = [68,42,132,157]
RX = [68,132,132,68,68]
RY = [42,42,157,157,42]
vmin = ShortFrequency[ShortFrequency > 0].min()
vmax = ShortFrequency[ShortFrequency > 0].max()
print("maxmin value", vmin, vmax)

vmin = LongFrequency[LongFrequency > 0].min()
vmax = LongFrequency[LongFrequency > 0].max()
print("long maxmin value", vmin, vmax)

ShortFrequency[ShortFrequency ==-9999] = np.nan
LongFrequency[LongFrequency==-9999] = np.nan

fig = plt.figure(figsize=(10,6))

spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])
# Short Rains

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_ylabel("Belg",fontdict={'fontname':'Times New Roman','fontsize':16})
# ax1.axis('off')
cmap = plt.get_cmap("coolwarm")
boundaries = [1.5,2.5, 3.5, 4.5, 5.5,6.5,7.5,8.5,9.5,10.5]
norm = mpl.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
cax = ax1.imshow(ShortFrequency, cmap=cmap, norm=norm)
cmap.set_under('white')
cbar = plt.colorbar(cax, ax=ax1,fraction=0.0355, pad=0.04)
cbar.ax.get_yaxis().set_ticks([])
cbar.ax.get_yaxis().set_ticks([2,3,4,5,6,7,8,9,10])
cbar.ax.get_yaxis().set_ticklabels([2,3,4,5,6,7,8,9,10])
cbar.ax.get_yaxis().set_tick_params(direction='in')
ax1.text(0.05, 0.95,'(a)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax1.transAxes)
ax1.plot(RX,RY,"r")
# Long Rains
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel("\nDrought Frequency",fontdict={'fontname':'Times New Roman','fontsize':16})
ax3.set_ylabel("Kiremit",fontdict={'fontname':'Times New Roman','fontsize':16})
cax2 = ax3.imshow(LongFrequency, cmap=cmap, norm=norm)
cmap.set_under('white')
cbar2 = plt.colorbar(cax2,ax=ax3,fraction=0.0355, pad=0.04)
cbar2.ax.get_yaxis().set_ticks([])
cbar2.ax.get_yaxis().set_ticks([2,3,4,5,6,7,8,9,10])
cbar2.ax.get_yaxis().set_ticklabels([2,3,4,5,6,7,8,9,10])
cbar2.ax.get_yaxis().set_tick_params(direction='in')
ax3.text(0.05, 0.95,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax3.transAxes)
ax3.plot(RX,RY,"r")

# **********


ShortFrequency,ShortPMat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage","Short")
LongFrequency,LongPMat = AnomalyFrequency(r"D:\Cornell\EthiopianDrought\0ExperimentData\Precipitation_Data\P_SeasonlyAverageImage","Long")

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
print(np.array(SNWTS).mean())
print(np.array(LNWTS).mean())
ax2.set_xticks([])
ax2.set_xlabel("P(mm)",fontdict={'fontname':'Times New Roman','fontsize':16})
ax2.plot(range(Years),SNWTS,'k-*',label="precipitation")
ax2.set_xticks(range(0,39,5),[str(1981+y)[2:] for y in range(0,39,5)])
ax2.text(0.05, 0.95,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax2.transAxes)


ax4.set_xlabel("Year(1981-2018)",fontdict={'fontname':'Times New Roman','fontsize':16})
ax4.set_xlabel("P(mm)",fontdict={'fontname':'Times New Roman','fontsize':16})
ax4.plot(range(Years),LNWTS,'k-*',label="precipitation")

ax4.set_xticks(range(0,39,5))
ax4.set_xticklabels([str(1981+y) for y in range(0,39,5)])
# plt.legend(loc=2)
ax4.text(0.05, 0.95,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax4.transAxes)




plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.05)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_1_Long.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()







