from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
# ET\
def ExtractRainFall(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\Chirps2"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])] for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []


    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(path,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        assert os.path.exists(chirps_file),"the file does not exists {}".format(chirps_file)
        data = gdal.Open(chirps_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values

def ExtractEVI(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        evi_file = os.path.join(path,
                                   str(dt.year) + "." + str(dt.month).zfill(2) + ".01"+".tif")
        assert os.path.exists(evi_file),"the file does not exists {}".format(evi_file)
        data = gdal.Open(evi_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values
def ExtractNEWSIF(CE,CN):
    path = r"D:\Cornell\NewSIF005Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_200208.nc.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path,"SIF005_{}{}.nc.tif".format(str(dt.year),str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file),"the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values
def ExtractGOSIF(CE,CN):
    path = r"D:\Cornell\GOSIFV002Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path,"GOSIF_{}.M{}.tif".format(str(dt.year),str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file),"the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values
def ExtractET(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\AData\ETClip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ETClip\ET.v3.3a198001.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path,"ET.v3.3a{}{}.tif".format(str(dt.year),str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file),"the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values
def ExtractSM(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5\ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-199101.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path,"ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-{}{}.tif".format(str(dt.year),str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file),"the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values
def ExtractPVI(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\AData\PVI"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI\pvi_2010.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = [[float(CE[i]),float(CN[i])]for i in range(len(CE))]
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]
    Col = np.array(Col)
    Row = np.array(Row)
    XSize ,YSize = raster.RasterXSize,raster.RasterYSize
    mask = (Col >= 0) & (Row >=0) & (Col < XSize) & (Row < YSize)
    assert any(mask) == True,"there is a - col and row"
    Col = Col.astype(np.int)
    Row = Row.astype(np.int)
    Ind = Row*XSize + Col
    Values = []

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path,"pvi_{}.tif".format(str(dt.year)))
        assert os.path.exists(newsif_file),"the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        Values.append(np.take(data,Ind))
    return Values

def ExtractCSV(path):
    CE = []
    CN = []
    CropList = []
    CropId = []
    CropType = os.path.basename(path)[:-4]
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            temp = []
            CE.append(row["CE"])
            CN.append(row["CN"])
            temp = []
            temp.append(row[CropType+"2010"])
            temp.append(row[CropType+"2011"])
            temp.append(row[CropType+"2012"])
            temp.append(row[CropType+"2013"])
            temp.append(row[CropType+"2014"])
            temp.append(row[CropType+"2015"])
            temp.append(row[CropType+"2016"])
            CropId.append(row["CropID"])
            CropList.append(temp)
    CropList = np.array(CropList)
    RFValues = ExtractRainFall(CE, CN)
    EVIValues = ExtractEVI(CE, CN)
    NEWSIFValues = ExtractNEWSIF(CE, CN)
    GOSIFValues = ExtractGOSIF(CE, CN)
    ETValues= ExtractET(CE, CN)
    SMValues = ExtractSM(CE, CN)
    PVIValues = ExtractPVI(CE, CN)

    Month =[1,2,3,4,5,6,7,8,9,10,11,12]
    Year =[2010,2011,2012,2013,2014,2015,2016]
    YM = []
    for year in Year:
        for m in Month:
            YM.append(str(year)+str(m).zfill(2))
    DataDict = {CropType + '2010': CropList[:, 0], CropType + '2011': CropList[:, 1], CropType + '2012': CropList[:, 2],
          CropType + '2013': CropList[:, 3], CropType + '2014': CropList[:, 4], CropType + '2015': CropList[:, 5],
          CropType + '2016': CropList[:, 6], "CE": CE, "CN": CN, "CropID": CropId}
    for idx, ym in enumerate(YM):
        DataDict["RF"+ym] = RFValues[idx]
    for idx, ym in enumerate(YM):
        DataDict["EVI"+ym] = EVIValues[idx]
    for idx, ym in enumerate(YM):
        DataDict["NEWSIF"+ym] = NEWSIFValues[idx]
    for idx, ym in enumerate(YM):
        DataDict["GOSIF"+ym] = GOSIFValues[idx]
    for idx, ym in enumerate(YM):
        DataDict["ET"+ym] = ETValues[idx]
    for idx, ym in enumerate(YM):
        DataDict["SM"+ym] = SMValues[idx]
    for idx, ym in enumerate(Year):
        DataDict["PVI"+str(ym)] = PVIValues[idx]
    print(DataDict)
    df = pd.DataFrame(DataDict)
    outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\{}_From2010-2016.csv".format(CropType)
    df.to_csv(outpath, index=False)


# croppath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\WHEATOPH.csv"
# ExtractCSV(croppath)

def LogOutRain(CE,CN,flag):
    path = r"D:\Cornell\EthiopianDrought\Chirps2"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    p0,p1 = ct.TransformPoint(float(CE),float(CN))[0:2]
    Col = int((p1-geot[0])/geot[1])
    Row = int((p0-geot[3])/geot[5])
    print(p1,p0)
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(path,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        assert os.path.exists(chirps_file),"the file does not exists {}".format(chirps_file)
        data = gdal.Open(chirps_file).ReadAsArray()

        print(dt.date(),Row,Col,data[Row][Col])
def LogOutEVI(CE,CN):
    path = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    geot = raster.GetGeoTransform()
    proj = osr.SpatialReference()

    proj.ImportFromWkt(str(raster.GetProjection()))
    print(proj)
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)

    p0,p1 = ct.TransformPoint(CE,CN)[0:2]
    Col = int((p1-geot[0])/geot[1])
    Row = int((p0-geot[3])/geot[5])
    print(p1,p0,Col,Row)

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        evi_file = os.path.join(path,
                                   str(dt.year) + "." + str(dt.month).zfill(2) + ".01"+".tif")
        assert os.path.exists(evi_file),"the file does not exists {}".format(evi_file)
        data = gdal.Open(evi_file).ReadAsArray()
        print(dt.date(),Row,Col,data[Row][Col])

LogOutRain(516715.929492471,1462445.57541472,"")
# LogOutEVI(663645.968351055,987489.22728797)
# 38.386258233803815 14.329207869071974

