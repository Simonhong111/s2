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
# path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask60.csv"
#

path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask60.csv"
#

DataType ="Mask"
Year = "2018"
fig = plt.figure(figsize=(10, 10))
# mMonthType = "Long"
plt.title(Year + " For {}\n".format(DataType), fontsize=16)


global point
point = []
def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    point.append([event.xdata,event.ydata])

titlename = "{} LongRains NSIF vs PVI".format(Year)
Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var3List)
xlabel = "NEWSIF"

scatter = plt.scatter(Var1List, Var3List)
titles = titlename
plt.title(titles, fontsize=10)
plt.xlabel(xlabel)
plt.ylabel("PVI")
plt.text(0.5, 0.4, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
plt.text(0.5, 0.35, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

from osgeo import ogr


ring = ogr.Geometry(ogr.wkbLinearRing)
for p in point:
    ring.AddPoint(p[0],p[1])
ring.AddPoint(point[0][0],point[0][1])
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)


ROW2 = ROW.tolist()
COL2 = COL.tolist()
valid = []
print("point",point)
for id,r in enumerate(ROW2):

    point2 = ogr.Geometry(ogr.wkbPoint)
    point2.AddPoint(Var1List[id],Var3List[id])
    # print("p2",point2)
    if poly.Contains(point2) :
        print("ddd")
        valid.append(id)

R = ROW[valid]
C=COL[valid]


data = gdal.Open(r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip60.tif").ReadAsArray()
plt.imshow(data)
plt.scatter(C,R)
plt.show()

plt.scatter(Var1List,Var3List)
plt.scatter(Var1List[valid],Var3List[valid],s=2)
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid],Var3List[valid])

print(r_value)
