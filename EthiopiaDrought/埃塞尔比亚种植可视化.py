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
    return Points

def GetInfoFromFeature(feature):
    Geom = feature.GetGeometryRef()
    Area = Geom.GetArea()

    Centroid = [Geom.Centroid().GetX(),Geom.Centroid().GetY()]
    fcount = Geom.GetGeometryCount()
    gname = Geom.GetGeometryName()
    print("geometrytype",gname)
    Points = []

    if gname == "MULTIPOLYGON":
        for gid in range(fcount):
            Points.extend(Poy2Points(Geom.GetGeometryRef(gid)))
    elif gname =="POLYGON":
        Points.extend(Poy2Points(Geom))

    # proj = feature.GetGeometryRef().GetSpatialReference()
    RK_Code = feature.GetField(12)
    EA_Name = feature.GetField(13)
    print(RK_Code)
    print(EA_Name)
    if "-" not in EA_Name:
        EA_Code = str(RK_Code).split(".")[0] + str(EA_Name)[-2:]
    else:
        EA_Code = str(RK_Code).split(".")[0]+str(EA_Name).split("-")[1]

    return Area,Points,EA_Code,Centroid
def LogInfoOfFeature(feature):

    FieldCount = feature.GetFieldCount()
    for i in range(FieldCount):
        fieldName = feature.GetFieldDefnRef(i).GetName()
        field =  feature.GetField(i)
        print(fieldName,field)
    print("FID",feature.GetFID())


#
path = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(path)
layer = dataset.GetLayer()
feature = layer.GetFeature(7760)
LogInfoOfFeature(feature)

Area,Points,EA_Code,Centroid = GetInfoFromFeature(feature)
print(Area,Area/10000)
Points = np.array(Points)

X,Y = Points[:,0],Points[:,1]
print(Centroid)
plt.plot(X,Y)
# plt.scatter(X[0:1],Y[0:1],c="r")
# plt.scatter(X[-1:],Y[-1:],c="g")
plt.scatter(Centroid[0],Centroid[1],c='r',s=200)
plt.show()

aggpath = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"

HeadLines = ["EA_code_merge","AREAH2010","REGIONCODE","ZONECODE","WOREDACODE","KEBELECODE",
"WHEATAREA2010","WHEATOUTPUT2010","MAIZEAREA2010","MAIZEOUTPUT2010",
"BARLEYAREA2010","BARLEYOUTPUT2010","SORGHUMAREA2010","SORGHUMOUTPUT2010",
"TEFFAREA2010","TEFFOUTPUT2010","HOUSEHOLD2010",
"AREAH2011",
"WHEATAREA2011","WHEATOUTPUT2011","MAIZEAREA2011","MAIZEOUTPUT2011",
"BARLEYAREA2011","BARLEYOUTPUT2011","SORGHUMAREA2011","SORGHUMOUTPUT2011",
"TEFFAREA2011","TEFFOUTPUT2011","HOUSEHOLD2011","AREAH2012",
"WHEATAREA2012","WHEATOUTPUT2012","MAIZEAREA2012","MAIZEOUTPUT2012",
"BARLEYAREA2012","BARLEYOUTPUT2012","SORGHUMAREA2012","SORGHUMOUTPUT2012",
"TEFFAREA2012","TEFFOUTPUT2012","HOUSEHOLD2012",
"AREAH2013",
"WHEATAREA2013","WHEATOUTPUT2013","MAIZEAREA2013","MAIZEOUTPUT2013",
"BARLEYAREA2013","BARLEYOUTPUT2013","SORGHUMAREA2013","SORGHUMOUTPUT2013",
"TEFFAREA2013","TEFFOUTPUT2013","HOUSEHOLD2013",
"AREAH2014",
"WHEATAREA2014","WHEATOUTPUT2014","MAIZEAREA2014","MAIZEOUTPUT2014",
"BARLEYAREA2014","BARLEYOUTPUT2014","SORGHUMAREA2014","SORGHUMOUTPUT2014",
"TEFFAREA2014","TEFFOUTPUT2014","HOUSEHOLD2014",
"AREAH2015",
"WHEATAREA2015","WHEATOUTPUT2015","MAIZEAREA2015","MAIZEOUTPUT2015",
"BARLEYAREA2015","BARLEYOUTPUT2015","SORGHUMAREA2015","SORGHUMOUTPUT2015",
"TEFFAREA2015","TEFFOUTPUT2015","HOUSEHOLD2015",
"AREAH2016",
"WHEATAREA2016","WHEATOUTPUT2016","MAIZEAREA2016","MAIZEOUTPUT2016",
"BARLEYAREA2016","BARLEYOUTPUT2016","SORGHUMAREA2016","SORGHUMOUTPUT2016",
"TEFFAREA2016","TEFFOUTPUT2016","HOUSEHOLD2016",
"WHEATOPH2010","WHEATOPH2011","WHEATOPH2012","WHEATOPH2013","WHEATOPH2014","WHEATOPH2015",
"WHEATOPH2016","MAIZEOPH2010","MAIZEOPH2011","MAIZEOPH2012","MAIZEOPH2013","MAIZEOPH2014",
"MAIZEOPH2015","MAIZEOPH2016","BARLEYOPH2010","BARLEYOPH2011","BARLEYOPH2012","BARLEYOPH2013",
"BARLEYOPH2014","BARLEYOPH2015","BARLEYOPH2016","SORGHUMOPH2010","SORGHUMOPH2011","SORGHUMOPH2012",
"SORGHUMOPH2013","SORGHUMOPH2014","SORGHUMOPH2015","SORGHUMOPH2016","TEFFOPH2010","TEFFOPH2011",
"TEFFOPH2012","TEFFOPH2013","TEFFOPH2014","TEFFOPH2015","TEFFOPH2016","_ID"]


def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
def CropYieldFromName(aggpath,CropType):
    CropType = [CropType + str(year) for year in range(2010, 2017)]
    ID = []
    Value = []
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
# feature = layer.GetFeature(35526)
# Area, Points, EA_Code, Centroid = GetInfoFromFeature(feature)
# print(Area, Points, EA_Code, Centroid)


CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
crtype = CTp[0]
CropList1,CropId1 = CropYieldFromName(aggpath,crtype)
CE = []
CN = []
CropList=[]
CropId=[]
for id,term in enumerate(CropId1):
    if int(term) == 58125:
        continue
    _ID = int(term)-1
    print("**",term)
    feature = layer.GetFeature(_ID)
    Area, Points, EA_Code, Centroid = GetInfoFromFeature(feature)
    assert term !=EA_Code,"the EA_Cord is wrong {}".format(term)
    CE.append(float(Centroid[0]))
    CN.append(float(Centroid[1]))
    CropList.append(CropList1[id])
    print(CropList1[id])
    CropId.append(term)

CropList = np.array(CropList)

# df = pd.DataFrame({crtype+'2010': CropList[:,0],crtype+'2011': CropList[:,1],crtype+'2012': CropList[:,2],
#                    crtype+'2013': CropList[:,3],crtype+'2014': CropList[:,4],crtype+'2015': CropList[:,5],
#                    crtype+'2016': CropList[:,6],"CE":CE,"CN":CN,"CropID":CropId})
# outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\{}.csv".format(crtype)
# df.to_csv(outpath,index=False)











