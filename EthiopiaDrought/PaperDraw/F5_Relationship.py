from osgeo import gdal,osr,ogr
import os
import numpy as np
from scipy import stats
from matplotlib import  pyplot as plt
import matplotlib.gridspec as gridspec
def CMIP5Comp(CMIP5Dir,SeasonType="Short"):

    YearNum = 12
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0
    for year in range(2007,2019):

        cm_file = os.path.join(CMIP5Dir,
                                   "cmip5_{}_{}.tif".format(SeasonType, str(year)))
        cm_data = gdal.Open(cm_file).ReadAsArray()

        Multidarr[:, :, band_id] = cm_data
        Mask[cm_data == -9999] = -9999
        band_id += 1

    return Multidarr,Mask

def ChirpsComp(CHIRPSDir,SeasonType="Short"):

    YearNum = 12
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0

    for year in range(2007,2019):

        ch_file = os.path.join(CHIRPSDir,
                                   "chirps-v2.0.{}_{}.tif".format(SeasonType, str(year)))

        ch_data = gdal.Open(ch_file).ReadAsArray()

        Multidarr[:, :, band_id] = ch_data
        Mask[ch_data == -9999] = -9999
        band_id += 1


    return Multidarr,Mask

def CMIP5PVIComp(CMPDir,SeasonType="Short"):

    YearNum = 12
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros((20, 13, YearNum))
    band_id = 0

    for year in range(2007,2019):

        cmp_file = os.path.join(CMPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        cmp_data = gdal.Open(cmp_file).ReadAsArray()

        Multidarr[:, :, band_id] = cmp_data
        Mask[cmp_data == -9999] = -9999
        band_id += 1

    return Multidarr,Mask

def ChirpsPVIComp(CHPDir,SeasonType="Short"):

    YearNum = 12
    Mask = np.zeros((20, 13))
    Multidarr = np.zeros(shape=(20, 13, YearNum),dtype=np.float)
    band_id = 0

    for year in range(2007,2019):

        chp_file = os.path.join(CHPDir,
                                   "{}_pvi_{}.tif".format(SeasonType, str(year)))

        chp_data = gdal.Open(chp_file).ReadAsArray()

        Multidarr[:, :, band_id] = chp_data
        Mask[chp_data == -9999] = -9999
        band_id += 1


    return Multidarr,Mask

def correlation(Multiarr1,Multiarr2):
    H,W,Chanel = Multiarr1.shape
    Coefficient = np.zeros(shape=(H,W),dtype=np.float)
    Pvalue = np.zeros(shape=(H,W),dtype=np.float)
    Flatten1 = Multiarr1.reshape(H*W,Chanel)
    Flatten2 = Multiarr2.reshape(H*W,Chanel)
    for row in range(H):
        for col in range(W):
            index = row*W + col
            slope, intercept, r_value, p_value, std_err = stats.linregress(Flatten1[index], Flatten2[index])
            Coefficient[row,col] = r_value
            Pvalue[row,col] = p_value

    return Coefficient,Pvalue



ref_path = r"D:\Cornell\EthiopianDrought\AData\CMIP5PVI\Big\long_pvi_2006.tif"
ref_raster = gdal.Open(ref_path)
geo_t = ref_raster.GetGeoTransform()

# 计算矢量边界

daShapefile = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()
geo = str(geo).split("((")[1].split("))")[0].split(",")
x = []
y = []
for term in geo:
    x.append(float(term.split(" ")[0]))
    y.append(float(term.split(" ")[1]))

x = np.array(x)
y = np.array(y)
x = (x - geo_t[0]) / geo_t[1]
y = (y - geo_t[3]) / geo_t[5]



yy = "2006-2018"
mm = 'Short'
cm_rf_path = r"D:\Cornell\EthiopianDrought\CMIPMonth\Big"
cm_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CMIP"
ch_rf_path = r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth\Big"
ch_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CHIRPS"


chrf,chrfMask = ChirpsComp(ch_rf_path,mm)
chpvi,chpviMask = ChirpsPVIComp(ch_pvi_path,mm)
cmrf,cmrfMask = CMIP5Comp(cm_rf_path,mm)
cmpvi,cmpviMask = CMIP5PVIComp(cm_pvi_path,mm)



rfCoef,rfPvalue = correlation(chrf,cmrf)
mask = (chrfMask == -9999) | (cmrfMask == -9999) | (rfPvalue <0.05)
print("short",rfCoef[mask == False].min(),rfCoef[mask == False].max())
rfCoef[mask] = np.nan

pviCoef,pviPvalue = correlation(chpvi,cmpvi)
mask2 = (chpviMask == -9999) | (cmpviMask == -9999) | (pviPvalue <0.05)
print("short pvi",pviCoef[mask2 == False].min(),pviCoef[mask2 == False].max())
pviCoef[mask2] = np.nan

fraction =0.069
fig = plt.figure(figsize=(10,6))
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])


cmap=plt.get_cmap("bwr")
vmin = -0.6
vmax = 0.6
ax1.text(0.95, 0.05,'(a)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax1.transAxes)

cax1 = ax1.imshow(rfCoef, cmap=cmap, vmin=vmin, vmax=vmax)
cbar1 = plt.colorbar(cax1, ax=ax1, fraction=fraction, pad=0.04)
ax1.plot(x,y)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_ylabel("Belg",fontdict={'fontname':'Times New Roman','fontsize':16})


# chiprs p

ax2.text(0.95, 0.05,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax2.transAxes)

chrf[chrf == -9999] = np.nan
cax2 = ax2.imshow(pviCoef, cmap=cmap, vmin=vmin, vmax=vmax)
cbar2 = plt.colorbar(cax2, ax=ax2, fraction=fraction, pad=0.04)
ax2.plot(x,y)
ax2.set_xticks([])
ax2.set_yticks([])

#********************
yy = "2006-2018"
mm = 'Long'
cm_rf_path = r"D:\Cornell\EthiopianDrought\CMIPMonth\Big"
cm_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CMIP"
ch_rf_path = r"D:\Cornell\EthiopianDrought\ChirpsDailyMonth\Big"
ch_pvi_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\PVI_Data\PVI_4M_5D_CHIRPS"

chrf,chrfMask = ChirpsComp(ch_rf_path,mm)
chpvi,chpviMask = ChirpsPVIComp(ch_pvi_path,mm)
cmrf,cmrfMask = CMIP5Comp(cm_rf_path,mm)
cmpvi,cmpviMask = CMIP5PVIComp(cm_pvi_path,mm)

rfCoef,rfPvalue = correlation(chrf,cmrf)
mask = (chrfMask == -9999) | (cmrfMask == -9999) | (rfPvalue <0.05)
print("long",rfCoef[mask == False].min(),rfCoef[mask == False].max())
rfCoef[mask] = np.nan

pviCoef,pviPvalue = correlation(chpvi,cmpvi)
mask2 = (chpviMask == -9999) | (cmpviMask == -9999) | (pviPvalue <0.05)
print("long pvi",pviCoef[mask2 == False].min(),pviCoef[mask == False].max())
pviCoef[mask2] = np.nan

ax3.text(0.95, 0.05,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax3.transAxes)

cax1 = ax3.imshow(rfCoef, cmap=cmap, vmin=vmin, vmax=vmax)
cbar1 = plt.colorbar(cax1, ax=ax3, fraction=fraction, pad=0.04)
ax3.plot(x,y)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_ylabel("Kiremit",fontdict={'fontname':'Times New Roman','fontsize':16})
ax3.set_xlabel("Cmip5 P",fontdict={'fontname':'Times New Roman','fontsize':16})
print("d",np.array(x).min(),np.array(x).max())

# chiprs p

ax4.text(0.95, 0.05,'(f)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='right',
     verticalalignment='center',
     transform = ax4.transAxes)


cax2 = ax4.imshow(pviCoef, cmap=cmap, vmin=vmin, vmax=vmax)

cbar2 = plt.colorbar(cax2, ax=ax4, fraction=fraction, pad=0.04)
ax4.plot(x,y)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlabel("Chirps P",fontdict={'fontname':'Times New Roman','fontsize':16})
# cmip5 pvi

plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_5.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()






