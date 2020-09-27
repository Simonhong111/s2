# import os
# import numpy as np
# from matplotlib import pyplot as plt
# from osgeo import gdal
# strlist = "114.911396628747 30.7162636630687 114.932223205884 29.7260490255611 116.06723105416 29.7389165956819 116.057829637596 30.7296502513598 114.911396628747 30.7162636630687"
# strlist = strlist.split()
# x = [float(strlist[2*i]) for i in range(int(len(strlist)/2))]
# y = [float(strlist[2*i+1]) for i in range(int(len(strlist)/2))]
# x = np.array(x)
# y = np.array(y)
# # x = (x - 199980)/60
# # y = (3400020 - y)/60
# raster = gdal.Open(r"D:\aria2\hongboshi\S2A_MSIL2A_20190418T030551_N0211_R075_T50RKU_20190418T073633.SAFE\GRANULE\L2A_T50RKU_A019949_20190418T030639\IMG_DATA\R60m\T50RKU_20190418T030551_TCI_60m.jp2").ReadAsArray().astype(np.float)
# img = np.zeros((1830,1830,3),dtype=np.int)
# print(raster.shape)
# img[:,:,0] = raster[0,:,:]
# img[:,:,1] = raster[1,:,:]
# img[:,:,2] = raster[2,:,:]
# print(img.max())
# ax = plt.gca()
# # ax.imshow(img)
# # ax.invert_yaxis()
# ax.plot(y,x)
#
# plt.show()


from osgeo import osr,gdal
import numpy as np
crs = osr.SpatialReference()
crs.ImportFromEPSG(32650)
wgs = osr.SpatialReference()
wgs.ImportFromEPSG(4326)
ct = osr.CoordinateTransformation(wgs,crs)

points = [[30.956,113.468],[30.2374,113.468],[30.2374,114.508],[30.956,114.516],[30.956,113.468]]
win = ct.TransformPoints(points)
print(win)
E = np.array([term[0] for term in win])
N = np.array([term[1] for term in win])
print("*",E)
from matplotlib import  pyplot as plt

raster = gdal.Open(r"C:\Users\zmhwh\Desktop\Temp\t49RGQ_20191106t025909_B02_20m_mosaic.tif")
gt = raster.GetGeoTransform()
# print(raster.GetProjection())
COL = (E-gt[0])/gt[1]
ROW = (N - gt[3])/gt[5]
print(gt)
print(COL)
plt.imshow(raster.ReadAsArray())
plt.plot(ROW,COL,'r')
plt.show()