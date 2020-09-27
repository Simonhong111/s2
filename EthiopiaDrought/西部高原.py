import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import  gdal,osr,ogr
import shapefile as shp
from osgeo import ogr


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
    # print(data[key1].to_numpy())
    # print(data[key1].to_numpy().mean(axis=1))
    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    ROW = data["ROW"].to_numpy()
    COL = data["COL"].to_numpy()

    return Var1List,Var2List,PVIList,ROW,COL

ref_path =r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip60.tif"
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

daShapefile = r"D:\Cornell\EthiopianDrought\cropland\CropLand.shp"
# daShapefile = r"D:\Cornell\EthiopianDrought\northhighlandshp\NorthHighland.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()
geo = str(geo).split("((")[1].split("))")[0].split(",")
x2 = []
y2 = []
for term in geo:
    x2.append(float(term.split(" ")[0]))
    y2.append(float(term.split(" ")[1]))
print(x2)
print(y2)
x2 = np.array(x2)
y2 = np.array(y2)
x2 = (x2 - geo_t[0]) / geo_t[1]
y2 = (y2 - geo_t[3]) / geo_t[5]
x2= x2.tolist()
y2=y2.tolist()

ring = ogr.Geometry(ogr.wkbLinearRing)
for id,p in enumerate(x2):
    ring.AddPoint(x2[id],y2[id])
ring.AddPoint(x2[0],y2[0])
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)




path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\10day\PolyGonAgg_Mask50.csv"
DataType ="Mask"
Year = "2012"
fig = plt.figure(figsize=(10, 10))
plt.title(Year + " For {}\n".format(DataType), fontsize=16)


titlename = "{} LongRains NSIF vs PVI".format(Year)
Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Short", "NSIF", "RF", Year)


ROW2 = ROW.tolist()
COL2 = COL.tolist()
valid = []

for id,r in enumerate(ROW2):

    point2 = ogr.Geometry(ogr.wkbPoint)
    point2.AddPoint(COL2[id],ROW2[id])
    # print("p2",point2)
    if poly.Contains(point2) :
        print("ddd")
        valid.append(id)

R = ROW[valid]
C=COL[valid]


# data = gdal.Open(r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip60.tif").ReadAsArray()
# plt.imshow(data)
# plt.scatter(C,R)
# plt.plot(x,y)
# plt.show()

# plt.scatter(Var1List,Var3List)
# plt.scatter(Var1List[valid],Var3List[valid],s=2)
# plt.show()
for Year in range(2003,2019):
    Var1List, Var2List, Var3List, ROW, COL = monthAnalysis(path, "Long", "NSIF", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])

    print(Year,r_value)
