import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import gdal,osr,ogr
path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\5day\PolyGonAgg_Mask50.csv"

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

    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    ROW = data["ROW"].to_numpy()
    COL = data["COL"].to_numpy()
    return Var1List, Var2List, PVIList, ROW, COL
# monthAnalysis(path,"Short","EVI","RF","2010")


# 生成有效范围
ref_path =r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip60.tif"
ref_raster = gdal.Open(ref_path)
geo_t = ref_raster.GetGeoTransform()
daShapefile = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropSHP\CropLandR.shp"
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




path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\5day\PolyGonAgg_Mask50.csv"
DataType ="Mask"
Year = "2012"
fig = plt.figure(figsize=(10, 10))
plt.title(Year + " For {}\n".format(DataType), fontsize=16)


titlename = "{} LongRains NSIF vs PVI".format(Year)
Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", Year)


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


DataType ="Crop Mask"
Years = [str(year) for year in range(2003,2019)]
for Year in Years:
    fig = plt.figure(figsize=(15, 10))
    mMonthType = "Long"
    plt.title(Year + " For {}\n".format(DataType), fontsize=16)
    plt.xticks([])
    plt.yticks([])
    for m in range(1, 7):

        ax = fig.add_subplot(2, 3, m)
        if m == 1:
            titlename = "{} ShortRains EVI vs RainFall".format(Year)
            Var1List, Var2List, Var3List ,ROW,COL= monthAnalysis(path, "Short", "EVI", "RF", Year)
            xlabel = "EVI"

            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
        if m == 2:
            titlename = "{} ShortRains GOSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List ,ROW,COL= monthAnalysis(path, "Short", "GSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
            xlabel = "GOSIF"
        if m == 3:
            titlename = "{} ShortRains NewSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Short", "NSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
            xlabel = "NEWSIF"
        if m == 4:
            titlename = "{} LongRains EVI vs RainFall".format(Year)
            Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Long", "EVI", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
            xlabel = "EVI"
        if m == 5:
            titlename = "{} LongRains GOSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List ,ROW,COL= monthAnalysis(path, "Long", "GSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
            xlabel = "GOSIF"
        if m == 6:
            titlename = "{} LongRains NSIF vs RainFall".format(Year)
            Var1List, Var2List, Var3List,ROW,COL = monthAnalysis(path, "Long", "NSIF", "RF", Year)
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
            xlabel = "NEWSIF"
        Var1List = Var1List[valid]
        Var2List = Var2List[valid]
        scatter = ax.scatter(Var1List, Var2List)
        titles = titlename
        ax.set_title(titles, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("RainFall")
        if m <= 3:
            ax.text(0.1, 190, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            # ax.text(0.1, 180, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")
        else:

            ax.text(0.5, 380, "R = {}".format(round(r_value, 3)), fontsize=16, color="r")
            ax.text(0.5, 360, "P = {}".format(round(p_value, 3)), fontsize=16, color="r")

    fig.tight_layout()  # 调整整体空白
    path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\5day\Image5west", DataType + Year + "NRF50.jpg")
    plt.savefig(path2)
    plt.close()
    # plt.show()

#
#
