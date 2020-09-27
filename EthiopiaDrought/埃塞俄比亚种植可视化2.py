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
def CropYieldFromName(aggpath,mCropType):
    CropType = [mCropType + str(year) for year in range(2010, 2017)]

    CropArea = [mCropType[:-3]+"AREA" + str(year) for year in range(2010, 2017)]
    ID = []
    Value = []
    Area = []
    with open(aggpath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ID.append(row["_ID"])
            temp = []
            temp.append(row[CropType[0]])
            temp.append(row[CropType[1]])
            temp.append(row[CropType[2]])
            temp.append(row[CropType[3]])
            temp.append(row[CropType[4]])
            temp.append(row[CropType[5]])
            temp.append(row[CropType[6]])
            Value.append(temp)

            Atemp = []
            Atemp.append(row[CropArea[0]])
            Atemp.append(row[CropArea[1]])
            Atemp.append(row[CropArea[2]])
            Atemp.append(row[CropArea[3]])
            Atemp.append(row[CropArea[4]])
            Atemp.append(row[CropArea[5]])
            Atemp.append(row[CropArea[6]])
            Area.append(Atemp)

    CropList = []
    AreaList = []
    CropId = []
    for id in range(len(Value)):
        # print(Value[id])
        flag = True
        for v in Value[id]:
            if is_number(v) == False:
                flag = False
        for v in [id]:
            if is_number(v) == False:
                flag = False
        if flag:
            CropList.append(Value[id])
            AreaList.append(Area[id])
            CropId.append(ID[id])
    return CropList,AreaList,CropId

path = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(path)
layer = dataset.GetLayer()
aggpath = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"
CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
crtype = CTp[5]
CropList1,AreaList1,CropId1 = CropYieldFromName(aggpath,crtype)
CE = []
CN = []
CropList=[]
AreaList = []
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
    AreaList.append(AreaList1[id])
    print(CropList1[id])
    CropId.append(term)

CropList = np.array(CropList)
AreaList = np.array(AreaList)
df = pd.DataFrame({crtype+'2010': CropList[:,0],crtype[:-3]+"AREA"+'2010': AreaList[:,0],
                   crtype+'2011': CropList[:,1],crtype[:-3]+"AREA"+'2011': AreaList[:,1],
                   crtype+'2012': CropList[:,2],crtype[:-3]+"AREA"+'2012': AreaList[:,2],
                   crtype+'2013': CropList[:,3],crtype[:-3]+"AREA"+'2013': AreaList[:,3],
                   crtype+'2014': CropList[:,4],crtype[:-3]+"AREA"+'2014': AreaList[:,4],
                   crtype+'2015': CropList[:,5],crtype[:-3]+"AREA"+'2015': AreaList[:,5],
                   crtype+'2016': CropList[:,6],crtype[:-3]+"AREA"+'2016': AreaList[:,6],
                   "CropID":CropId})
outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\{}_AREA.csv".format(crtype)
df.to_csv(outpath,index=False)














