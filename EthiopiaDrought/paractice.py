import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import  gdal,osr,ogr
import shapefile as shp


def monthAnalysis(path,monthtype,var1,var2,year):
    if monthtype == "Short":
        months = [2,3,4,5]
    if monthtype == "Long":
        months =[6,7,8,9]
    data = pd.read_csv(path)
    key1 = [var1+str(year)+str(months[0]).zfill(2),var1+str(year)+str(months[1]).zfill(2),
            var1+str(year)+str(months[2]).zfill(2),var1+str(year)+str(months[3]).zfill(2)]
    key2 = [var2 + str(year) + str(months[0]).zfill(2), var2 + str(year) + str(months[1]).zfill(2),
            var2 + str(year) + str(months[2]).zfill(2), var2 + str(year) + str(months[3]).zfill(2)]
    print(data[key1].to_numpy())
    print(data[key1].to_numpy().mean(axis=1))
    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    ROW = data["ROW"].to_numpy()
    COL = data["COL"].to_numpy()

    return Var1List,Var2List,PVIList,ROW,COL
path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask60.csv"
#

# DataType ="Mask"
# Year = "2008"
# fig = plt.figure(figsize=(10, 10))
# # mMonthType = "Long"
# plt.title(Year + " For {}\n".format(DataType), fontsize=16)
# plt.xticks([])
# plt.yticks([])
#
#
# titlename = "{} LongRains NSIF vs PVI".format(Year)
# Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", Year)
# slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var3List)
# xlabel = "NEWSIF"
#
# scatter = plt.scatter(Var1List, Var3List)
# titles = titlename
# plt.title(titles, fontsize=10)
# plt.xlabel(xlabel)
# plt.ylabel("PVI")
# plt.text(0.5, 0.4, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
# plt.text(0.5, 0.35, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
# plt.show()



from osgeo import ogr

# Create ring
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(0.50,0.32)
ring.AddPoint(0.27,0.133)
ring.AddPoint(0.39,0.011)
ring.AddPoint(0.65,0.086)
ring.AddPoint(0.50,0.32)
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

# Var1List,Var2List,PVIList,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", "2009")
# ROW2 = ROW.tolist()
# COL2 = COL.tolist()
# valid = []
# for id,r in enumerate(ROW):
#
#     p1x,p1y = geo[0] + COL2[id] * geo[1], geo[3] + ROW2[id]*geo[5]
#     p2x, p2y = p1x + geo[1], p1y
#     p3x, p3y = p2x,p2y + geo[5]
#     p4x, p4y = p1x,p3y
#
#     ring = ogr.Geometry(ogr.wkbLinearRing)
#     ring.AddPoint(p1x,p1y)
#     ring.AddPoint(p2x,p2y)
#     ring.AddPoint(p3x,p3y)
#     ring.AddPoint(p4x,p4y)
#     ring.AddPoint(p1x,p1y)
#     poly = ogr.Geometry(ogr.wkbPolygon)
#     poly.AddGeometry(ring)
#
#
#     if not (geometry11.Intersect(poly) or geometry13.Intersect(poly)):
#         valid.append(id)


R = ROW[valid]
C=COL[valid]
data = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018\long_pvi_2008.tif").ReadAsArray()
plt.imshow(data)
plt.scatter(C,R)
plt.show()
# from shapely import wkt
# rainpath =r"D:\Cornell\EthiopianDrought\eth_rainpat\eth_rainpat.shp"
# driver = ogr.GetDriverByName('ESRI Shapefile')
# dataSource = driver.Open(rainpath, 0)
# layer = dataSource.GetLayer()
# feature1 = layer.GetFeature(1)
# feature11 = layer.GetFeature(11)
# feature13 = layer.GetFeature(13)
# geometry1 = feature1.GetGeometryRef()
# geometry11 = feature11.GetGeometryRef()
# geometry13 = feature13.GetGeometryRef()
#
#
# feature = layer.GetFeature(16)
# geometry = feature.GetGeometryRef()
# # 读取crop layer tif 影像
# crop = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVIDaily\long_pvi_2008.tif")
# geo = crop.GetGeoTransform()
# crop_prj = osr.SpatialReference()
# crop_prj.ImportFromWkt(str(crop.GetProjection()))
# rainprj = osr.SpatialReference()
# rainprj.ImportFromEPSG(20137)
# geometry1.TransformTo(crop_prj)
# geometry1.SwapXY()
# geometry11.TransformTo(crop_prj)
# geometry11.SwapXY()
# geometry13.TransformTo(crop_prj)
# geometry13.SwapXY()
#
#
# geometry.TransformTo(crop_prj)
# geometry.SwapXY()
#
#
# Var1List,Var2List,PVIList,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", "2009")
# ROW2 = ROW.tolist()
# COL2 = COL.tolist()
# valid = []
# for id,r in enumerate(ROW):
#
#     p1x,p1y = geo[0] + COL2[id] * geo[1], geo[3] + ROW2[id]*geo[5]
#     p2x, p2y = p1x + geo[1], p1y
#     p3x, p3y = p2x,p2y + geo[5]
#     p4x, p4y = p1x,p3y
#
#     ring = ogr.Geometry(ogr.wkbLinearRing)
#     ring.AddPoint(p1x,p1y)
#     ring.AddPoint(p2x,p2y)
#     ring.AddPoint(p3x,p3y)
#     ring.AddPoint(p4x,p4y)
#     ring.AddPoint(p1x,p1y)
#     poly = ogr.Geometry(ogr.wkbPolygon)
#     poly.AddGeometry(ring)
#
#
#     if not (geometry11.Intersect(poly) or geometry13.Intersect(poly)):
#         valid.append(id)
# x = Var1List[valid]
# y = PVIList[valid]
# plt.imshow(crop.ReadAsArray())
# plt.scatter(COL[valid],ROW[valid])
# # plt.scatter(COL[valid][x<0.2],ROW[valid][x<0.2])
# plt.show()
# # plt.scatter(Var1List[valid],PVIList[valid],s=1)
# #
# # slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
# # print(r_value)
# # plt.show()
#
# #
#
# RV = []
# for year in range(2003,2019):
#     Var1List, Var2List, PVIList, ROW, COL = monthAnalysis(path, "Long", "NSIF", "RF", str(year))
#     x = Var1List[valid]
#     y = PVIList[valid]
#     slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#     RV.append(r_value*(-1))
# plt.grid(True)
# plt.plot(range(2003,2019),RV)
# plt.show()
# x = Var1List[valid]
# y = PVIList[valid]
# # plt.scatter(Var1List,PVIList)
# plt.scatter(Var1List[valid],PVIList[valid],s=1)
#
# slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
# print(r_value)
# plt.show()
#
#
#
#
#
#



