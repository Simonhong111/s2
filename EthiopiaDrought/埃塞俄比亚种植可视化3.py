from osgeo import gdal,osr,ogr
import numpy as np
import os
from matplotlib import pyplot as plt
import csv
import pandas as pd
def Poy2Points(Poly):
    Polystr = str(Poly)
    Polystr = Polystr.split("((")[1]
    Polystr = Polystr.split("))")[0]
    Polystr = Polystr.split(",")
    Points = []
    for p in Polystr:
        if ")" in p:
            p = p[:-1]
        if "(" in p:
            p = p[1:]
        Points.append([float(p.split(" ")[0]),float(p.split(" ")[1])])
    print("points",Points)
    return Points
def GetROWCoLFromFeature(polygon,geotransform):
    Points = Poy2Points(polygon)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(32637)
    des_srs = osr.SpatialReference()
    des_srs.ImportFromEPSG(4326)
    ct = osr.CoordinateTransformation(src_srs,des_srs)
    LonLat = ct.TransformPoints(Points)
    Lon = []
    Lat = []
    for LL in LonLat:
        Lat.append(LL[0])
        Lon.append(LL[1])

def GetInfoFromFeature(feature):
    Geom = feature.GetGeometryRef()
    fcount = Geom.GetGeometryCount()
    gname = Geom.GetGeometryName()
    print("geometrytype",gname,"fcount",fcount)
    Points = []

    if gname == "MULTIPOLYGON":
        for gid in range(fcount):
            Points.extend(Poy2Points(Geom.GetGeometryRef(gid)))
    elif gname =="POLYGON":
        Points.extend(Poy2Points(Geom))

    RK_Code = feature.GetField(12)
    EA_Name = feature.GetField(13)
    print(RK_Code)
    print(EA_Name)
    if "-" not in EA_Name:
        EA_Code = str(RK_Code).split(".")[0] + str(EA_Name)[-2:]
    else:
        EA_Code = str(RK_Code).split(".")[0]+str(EA_Name).split("-")[1]

    return Points,EA_Code
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
def CropYieldFromName(aggpath):
    CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]
    CropType = []
    AreaType = []
    for mctp in CTp:
        for year in range(2010,2017):
            CropType.append(mctp + str(year))
    Area = ["WHEATAREA", "MAIZEAREA", "BARLEYAREA", "SORGHUMAREA", "TEFFAREA"]
    for area in Area:
        for year in range(2010,2017):
            AreaType.append(area+str(year))

    ID = []
    Value = []
    with open(aggpath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ID.append(row["_ID"])
            temp = []
            areatemp = []
            for i in range(len(CropType)):
                temp.append(row[CropType[i]])
                temp.append(row[AreaType[i]])
            Value.append(temp)

    CropList = []
    CropId = []
    for id in range(len(Value)):
        # print(Value[id])
        flag = True
        for v in Value[id]:
            if is_number(v) == False:
                flag = False
        if flag:
            CropList.append(Value[id])
            CropId.append(ID[id])
    return CropList,CropId

path = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(path)
layer = dataset.GetLayer()
aggpath = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"



CropList1,CropId1 = CropYieldFromName(aggpath)
CropList=[]
CropId=[]
for id,term in enumerate(CropId1):
    if int(term) == 58125:
        continue
    _ID = int(term)-1
    print("**",term)
    feature = layer.GetFeature(_ID)
    Points, EA_Code = GetInfoFromFeature(feature)
    assert term !=EA_Code,"the EA_Cord is wrong {}".format(term)
    CropList.append(CropList1[id])
    print(CropList1[id])
    CropId.append(term)

CropList = np.array(CropList)

DataDict = {}
CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]
CropType = []
for mctp in CTp:
    for year in range(2010,2017):
        CropType.append(mctp + str(year))
        CropType.append(mctp[:-3]+"AREA"+str(year))


Length = len(CropType)
for i in range(Length):
    DataDict[CropType[i]] = CropList[:,i]

DataDict["ID"] = CropId



# df = pd.DataFrame({crtype+'2010': CropList[:,0],crtype+'2011': CropList[:,1],crtype+'2012': CropList[:,2],
#                    crtype+'2013': CropList[:,3],crtype+'2014': CropList[:,4],crtype+'2015': CropList[:,5],
#                    crtype+'2016': CropList[:,6],"CE":CE,"CN":CN,"CropID":CropId})

df = pd.DataFrame(DataDict)
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Crop_2010-2016.csv"
df.to_csv(outpath,index=False)














