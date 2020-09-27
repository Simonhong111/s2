# from osgeo import gdal,osr,ogr
# import sys,os
# from Sentinel2Reader import SenL2AReader
# import numpy as np
# from GridGen import *
# from Sen2Clip import *
# from geo2mapxy import *
# import time
# def dms2dec(degree, minute, second, sphere="N"):
#     m = minute / 60.0
#     s = second / 3600.0
#     return degree + m + s
#
#
#
# ppath = r"D:\S2A_MSIL2A_20190417T033541_N0211_R061_T48SWC_20190417T073942.SAFE\GRANULE\L2A_T48SWC_A019935_20190417T033835\IMG_DATA\R20m\s2clsscom.tif"
# ppath =r"D:\S2A_MSIL2A_20190417T033541_N0211_R061_T48SWC_20190417T073942.SAFE\GRANULE\L2A_T48SWC_A019935_20190417T033835\IMG_DATA\R10m\T48SWC_20190417T033541_B02_10m.jp2"
#
#
#
# scope =(105.16996666666667, 34.200719444444445,105.20851111111111, 34.169491666666666)
#
# #
# # start = time.time()
# # prosrs = osr.SpatialReference()
# # # prosrs.ImportFromWkt(o.GetProjection())
# # prosrs.ImportFromEPSG(32648)
# # geosrs = prosrs.CloneGeogCS()
# # # print(prosrs)
# # ct = osr.CoordinateTransformation(geosrs, prosrs)
# # for i in range(100):
# #
# #     coords = ct.TransformPoint(scope[0], scope[1])
# #     # print(coords)
# #     # print(coords[:2])
# # end = time.time()
# # print(end-start)
#
#
# import s2sphere as s2
#
# o = gdal.Open(ppath)
# print("**",o.GetGeoTransform())
# GT = o.GetGeoTransform()
#
# geoUpxy = np.array([GT[0],GT[3]])
# geoBxy =geoUpxy  + np.array([10*16384,-10*16384])
#
#
# geoxy1=np.array([499980,3800040])
# geoxy2 =np.array([609775,3690245])
#
# geost1x = (geoxy1[0] - geoUpxy[0]+1e-8)/(10*16384)
# geost1y = (geoUpxy[1]-geoxy1[1]+1e-8)/(10*16384)
# geost1 =np.array([geost1x,geost1y])
#
#
# geost2x = (geoxy2[0] - geoUpxy[0])/(10*16384)
# geost2y = (geoUpxy[1]-geoxy2[1])/(10*16384)
# geost2 =np.array([geost2x,geost2y])

# print(geost1)
# print(geost2)
#
# p1 = [s2.CellId.st_to_uv(geost1[1]),s2.CellId.st_to_uv(geost1[0])]
# p2 = [s2.CellId.st_to_uv(geost2[1]),s2.CellId.st_to_uv(geost2[0])]
# p1 = s2.face_uv_to_xyz(1,p1[0],p1[1])
# p2 = s2.face_uv_to_xyz(1,p2[0],p2[1])
# pp1 = s2.LatLng.from_point(p1)
# pp2 = s2.LatLng.from_point(p2)
# # point1 = s2.CellId.from_face_ij(1,p1[0],p1[1]).to_lat_lng()
# # point2 = s2.CellId.from_face_ij(1,p2[0],p2[1]).to_lat_lng()
# print(pp1)
# print(pp2)

# p1 = s2.LatLng.from_degrees(33, -122)
# p2 = s2.LatLng.from_degrees(33.1, -122.1)
#
# r = s2.RegionCoverer()
# r.max_level=14
# r.min_level=1
# start = time.time()
# cell_ids = r.get_covering(s2.LatLngRect.from_point_pair(pp1, pp2))
# for cell in cell_ids:
#     print(cell.face(),cell.level())
#     print(np.floor(s2.CellId.uv_to_st(cell.get_center_uv()[0])*16384),np.floor(s2.CellId.uv_to_st(cell.get_center_uv()[1])*16384))
# end = time.time()
# print(end-start)
# print(len(cell_ids))







# help(s2sphere)
# msphere = s2sphere.CellId
# print(msphere)
# print(msphere.st_to_ij(0.1))
# print(s2sphere.CellId.st_to_ij(0.1))
# start = time.time()
# r = s2sphere.RegionCoverer()
# p1 = s2sphere.LatLng.from_degrees(33, -122)
# p2 = s2sphere.LatLng.from_degrees(33.1, -122.1)
# cell_ids = r.get_covering(s2sphere.LatLngRect.from_point_pair(p1, p2))
# print(cell_ids)
# end = time.time()
# print(end-start)
# IJ = s2.CellId()
# IJ.MAX_LEVEL = 14
#
# IJ = IJ.from_face_ij(0,geost1[0],geost1[1])
# print(IJ.id())

import matplotlib.pyplot as plt
import numpy as np

cm = plt.cm.get_cmap('RdYlBu')
xy = range(100)
z = range(100)
sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=100, s=35, cmap=cm)
# sc =plt.scatter([1,2],[3,4],c=1, vmin=0, vmax=292, s=35, cmap=cm)
plt.colorbar(sc)
plt.show()


