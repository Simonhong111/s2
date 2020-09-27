from epdatapro20082015 import *
from matplotlib import pyplot as plt
import os
import numpy as np
import glob
from matplotlib import cm

months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

yy = "2017"
rp = r"D:\Cornell\EthiopianDrought\Chirps2"
rpp = r"D:\Cornell\EthiopianDrought\AData\Chirps2Pars"
ep = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
epp = r"D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyPars"
pvip = r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018"
pvipp = r"D:\Cornell\EthiopianDrought\AData\PviParms"
fig = plt.figure(figsize=(8, 8))
plt.title("Crop land Map" + '\n', fontsize=14)
#
# for m in range(1, 13):
#     path = os.path.join(rp,"chirps-v2.0.{}.{}.tif".format(yy,str(m).zfill(2)))
#     ramn = gdal.Open(path).ReadAsArray()
#     ax = fig.add_subplot(3, 4, m )
#     ax.set_title("Precipitaion_" + months[m - 1])
#     mask = np.where(ramn > -9999)
#     ramn[ramn == -9999] = np.nan
#     vmin = ramn[mask].min()
#     vmax = ramn[mask].max()
#     print("maxmin value", vmin, vmax)
#     cmap = plt.get_cmap("RdBu")
#
#     cax = ax.imshow(ramn, cmap=cmap, vmin=0, vmax=400)
#     cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
#     ax.set_xticks([])
#     ax.set_yticks([])


# for m in range(1, 13):
#     path = os.path.join(ep,"{}.{}.01.tif".format(yy,str(m).zfill(2)))
#
#     eamn = gdal.Open(path).ReadAsArray()*1.0
#     ax = fig.add_subplot(3, 4, m)
#     ax.set_title("EVI_" + months[m - 1])
#     mask = np.where(eamn > -3000)
#     eamn[eamn == -3000] = np.nan
#     eamn[mask] = eamn[mask]*0.0001
#     vmin = eamn[mask].min()
#     vmax = eamn[mask].max()
#     print("maxmin value",vmin,vmax)
#     cax = ax.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-0.3, vmax=0.9)
#     cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
#     ax.set_xticks([])
#     ax.set_yticks([])

# path = os.path.join(pvip,"pvi_{}.tif".format(yy))
# eamn = gdal.Open(path).ReadAsArray()
# ax = fig.add_subplot(1,3,1)
# ax.set_title("PVI_Yearly")
# mask = np.where(eamn > -9999)
# eamn[eamn == -9999] = np.nan
# vmin = eamn[mask].min()
# vmax = eamn[mask].max()
# print("maxmin value",vmin,vmax)
# cax = ax.imshow(eamn, cmap=plt.get_cmap("bwr"), vmin=0, vmax=0.6)
# cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
# ax.set_xticks([])
# ax.set_yticks([])
#
# path = os.path.join(pvip,"short_pvi_{}.tif".format(yy))
# ax = fig.add_subplot(1,3,2)
# eamn = gdal.Open(path).ReadAsArray()
# ax.set_title("PVI_Short_rains")
# mask = np.where(eamn > -9999)
# eamn[eamn == -9999] = np.nan
# vmin = eamn[mask].min()
# vmax = eamn[mask].max()
# print("maxmin value",vmin,vmax)
# cax = ax.imshow(eamn, cmap=plt.get_cmap("bwr"), vmin=0, vmax=0.6)
# cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
# ax.set_xticks([])
# ax.set_yticks([])
#
# path = os.path.join(pvip,"pvi_{}.tif".format(yy))
# eamn = gdal.Open(path).ReadAsArray()
# ax = fig.add_subplot(1,3,3)
# ax.set_title("PVI_Long_rains")
# mask = np.where(eamn > -9999)
# eamn[eamn == -9999] = np.nan
# vmin = eamn[mask].min()
# vmax = eamn[mask].max()
# print("maxmin value",vmin,vmax)
# cax = ax.imshow(eamn, cmap=plt.get_cmap("bwr"), vmin=0, vmax=0.6)
# cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
# ax.set_xticks([])
# ax.set_yticks([])
rain = r'D:\Cornell\MCD12C1V006Clip\MCD12C1.A2003001.006.2018053185458.hdf.tif'
rain = gdal.Open(rain)
geo_t = rain.GetGeoTransform()
print(geo_t)

daShapefile = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()

# for feature in layer:
#     geom = feature.GetGeometryRef()
#     print (geom,"*********************")

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

for id,year in enumerate([2008,2015,2016,2017]):
    path = glob.glob(os.path.join(r"D:\Cornell\MCD12C1V006Clip","MCD12C1.A{}001.006.*.hdf.tif".format(year)))[0]
    ramn = gdal.Open(path).ReadAsArray()*1.0
    ax = fig.add_subplot(2, 2, id+1)
    ax.set_title("{} cropland mask ".format(year))

    ramn[ramn != 12] = np.nan
    ax.imshow(ramn)
    ax.plot(x, y)
    ax.set_xticks([])
    ax.set_yticks([])
fig.tight_layout()  # 调整整体空白
plt.show()