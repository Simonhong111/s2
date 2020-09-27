from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt
# ET\

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
    # print("points",Points)
    return Points

def ExtractRainFall(feature,defval=-9999):
    path = r"D:\Cornell\EthiopianDrought\Chirps2"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop =  datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    invct = osr.CoordinateTransformation(proj,src)
    Centroid = [feature.Centroid().GetX(), feature.Centroid().GetY()]
    centroidY,centroidX = ct.TransformPoint(Centroid[0],Centroid[1])[0:2]
    Points = Poy2Points(feature)
    ccol = int((centroidX - geot[0]) / geot[1])
    crow = int((centroidY - geot[3]) / geot[5])
    ccpixel = ogr.Geometry(ogr.wkbLinearRing)
    ccpc1, ccpr1 = geot[0] + ccol * geot[1], geot[3] + crow * geot[5]
    ccpc2, ccpr2 = geot[0] + (ccol + 1) * geot[1], geot[3] + crow * geot[5]
    ccpc3, ccpr3 = geot[0] + (ccol + 1) * geot[1], geot[3] + (crow + 1) * geot[5]
    ccpc4, ccpr4 = geot[0] + ccol * geot[1], geot[3] + (crow + 1) * geot[5]

    cpjc1, cpjr1 = invct.TransformPoint(ccpr1, ccpc1)[0:2]
    cpjc2, cpjr2 = invct.TransformPoint(ccpr2, ccpc2)[0:2]
    cpjc3, cpjr3 = invct.TransformPoint(ccpr3, ccpc3)[0:2]
    cpjc4, cpjr4 = invct.TransformPoint(ccpr4, ccpc4)[0:2]

    ccpixel.AddPoint(cpjc1, cpjr1)
    ccpixel.AddPoint(cpjc2, cpjr2)
    ccpixel.AddPoint(cpjc3, cpjr3)
    ccpixel.AddPoint(cpjc4, cpjr4)
    ccpixel.AddPoint(cpjc1, cpjr1)
    ccpixelpoly = ogr.Geometry(ogr.wkbPolygon)
    ccpixelpoly.AddGeometry(ccpixel)

    if feature.Area() < 0.20*ccpixelpoly.Area():
        # print("are", feature.Area() / 1000000, ccpixelpoly.Area() / 1000000)
        return None,None,None,None,None,None,None,None

    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]


    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))


    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pc1,pr1 =  geot[0] + icol*geot[1],geot[3] + irow*geot[5]
            pc2, pr2 = geot[0] + (icol+1) * geot[1], geot[3] + irow * geot[5]
            pc3, pr3 = geot[0] + (icol+1) * geot[1], geot[3] + (irow+1) * geot[5]
            pc4, pr4 = geot[0] + icol * geot[1], geot[3] + (irow+1) * geot[5]

            pjc1,pjr1 = invct.TransformPoint(pr1,pc1)[0:2]
            pjc2, pjr2 = invct.TransformPoint(pr2,pc2)[0:2]
            pjc3, pjr3 = invct.TransformPoint(pr3,pc3)[0:2]
            pjc4, pjr4 = invct.TransformPoint(pr4,pc4)[0:2]

            pixel.AddPoint(pjc1,pjr1)
            pixel.AddPoint(pjc2,pjr2)
            pixel.AddPoint(pjc3,pjr3)
            pixel.AddPoint(pjc4,pjr4)
            pixel.AddPoint(pjc1,pjr1)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if feature.Intersect(pixelpoly):
                intersectionpoly = feature.Intersection(pixelpoly)
                if intersectionpoly.Area() >= ccpixelpoly.Area()*0.2 :

                    vCol.append(icol)
                    vRow.append(irow)

    # data = raster.ReadAsArray()
    if len(vCol) ==0:
        return None,None,None,None,None,None,None,None
    # print(vCol,vRow)
    ind = np.array(vRow)*XSize + np.array(vCol)


    EVI = []
    RF = []
    NSIF = []
    GSIF = []
    COLSet = []
    ROWSet = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        rf_file = os.path.join(path,"chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        gsif_file = os.path.join(r"D:\Cornell\GOSIFV002Clip","GOSIF_{}.M{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        nsif_file = os.path.join(r"D:\Cornell\NewSIF005Clip","SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))
        evi_file = os.path.join(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia","{}.{}.01.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(rf_file),"the file does not exists {}".format(rf_file)
        assert os.path.exists(gsif_file), "the file does not exists {}".format(gsif_file)
        assert os.path.exists(nsif_file), "the file does not exists {}".format(nsif_file)
        assert os.path.exists(evi_file), "the file does not exists {}".format(evi_file)
        # print(chirps_file)
        RFData = gdal.Open(rf_file).ReadAsArray()
        GSIFData = gdal.Open(gsif_file).ReadAsArray()
        NSIFData = gdal.Open(nsif_file).ReadAsArray()
        EVIData = gdal.Open(evi_file).ReadAsArray()


        RFV = np.take(RFData,ind)
        GSIFV = np.take(GSIFData, ind)
        NSIFV = np.take(NSIFData, ind)
        EVIV = np.take(EVIData, ind)

        EVI.append(EVIV.tolist())
        RF.append(RFV.tolist())
        NSIF.append(NSIFV.tolist())
        GSIF.append(GSIFV.tolist())
        del RFData
        del GSIFData
        del NSIFData
        del EVIData
        del RFV
        del GSIFV
        del NSIFV
        del EVIV
    COLSet.append(vCol)
    ROWSet.append(vRow)

    SPVI = []
    LPVI = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        short_file = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVI",
                                  "short_pvi_{}.tif".format(dt.year))
        long_file = os.path.join(r"D:\Cornell\EthiopianDrought\AData\PVI",
                                 "long_pvi_{}.tif".format(dt.year))

        assert os.path.exists(short_file), "the file does not exists {}".format(short_file)
        assert os.path.exists(long_file), "the file does not exists {}".format(long_file)


        ShortData = gdal.Open(short_file).ReadAsArray()
        LongData = gdal.Open(long_file).ReadAsArray()


        ShortV = np.take(ShortData, ind)
        LongV = np.take(LongData, ind)


        SPVI.append(ShortV.tolist())
        LPVI.append(LongV.tolist())

        del ShortData
        del LongData


    return EVI,RF,NSIF,GSIF,SPVI,LPVI,COLSet,ROWSet







def ExtractCSV(path):

    CropId = []
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            CropId.append(row["_ID"])

    shppath = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shppath)
    layer = dataset.GetLayer()
    EVID,RFD,NSIFD,GSIFD,SPVID,LPVID,COLSet,ROWSet= [],[],[],[],[],[],[],[]

    MYID = []
    for index, ID in enumerate(CropId):
        id = int(ID) - 1
        if id ==58124:
            continue
        feature = layer.GetFeature(id)
        Geom = feature.GetGeometryRef()
        print(id)
        fcount = Geom.GetGeometryCount()
        gname = Geom.GetGeometryName()

        if gname == "MULTIPOLYGON":
            continue
        elif gname == "POLYGON" and fcount == 1:

            EVI,RF,NSIF,GSIF,SPVI,LPVI,COL,ROW = ExtractRainFall(Geom)

            if EVI == None:
                continue
            else:
                MYID.append(ID)

    watershape_path = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"
    driver2 = ogr.GetDriverByName("ESRI Shapefile")
    dataset2 = driver2.Open(watershape_path)
    layer2 = dataset2.GetLayer()

    feature2 = layer2.GetFeature(0)
    geometry2 = feature2.GetGeometryRef()
    mPoint = Poy2Points(geometry2)
    mPoint = np.array(mPoint)
    print(mPoint)

    modis = gdal.Open(r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip.tif")

    modis_geo = modis.GetGeoTransform()
    modis_projection = modis.GetProjection()  # wgs84
    watershp_proj = osr.SpatialReference()
    watershp_proj.ImportFromEPSG(32637)
    modis_proj = osr.SpatialReference()
    modis_proj.ImportFromWkt(str(modis_projection))

    x = (mPoint[:, 0] -modis_geo[0])/modis_geo[1]
    y = (mPoint[:, 1] - modis_geo[3])/modis_geo[5]
    modis_raster = modis.ReadAsArray()*0
    modis_raster[modis_raster ==0] = np.nan
    plt.imshow(modis_raster)
    plt.plot(x,y)
    # plt.show()

    for index, ID in enumerate(MYID):
        id = int(ID) - 1
        feature = layer.GetFeature(id)
        Geom = feature.GetGeometryRef()
        Geom.AssignSpatialReference(watershp_proj)
        # Geom.SwapXY()

        Geom.TransformTo(modis_proj)
        # print("hh",Geom)

        mp = Poy2Points(Geom)
        mp = np.array(mp)

        x1 = (mp[:, 1] - modis_geo[0]) / modis_geo[1]
        y1 = (mp[:, 0] - modis_geo[3]) / modis_geo[5]
        plt.plot(x1,y1)
    plt.show()





croppath = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"


ExtractCSV(croppath)



