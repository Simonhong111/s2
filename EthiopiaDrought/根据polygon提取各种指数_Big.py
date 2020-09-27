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
        print("are", feature.Area() / 1000000, ccpixelpoly.Area() / 1000000)
        return None

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
                    weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()
    if len(vCol) ==0:
        return None
    # print(vCol,vRow)
    ind = np.array(vRow)*XSize + np.array(vCol)
    weigth = np.array(weigth)
    # pvalues = np.take(data,ind)

    # print(pvalues)
    # for id,_ in enumerate(vRow):
    #     print()
    #     print(data[vRow[id]][vCol[id]])
    # weigth = np.array(weigth)
    # vmask = np.where(pvalues != defval)
    # print(weigth[vmask],pvalues[vmask],weigth[vmask]*pvalues[vmask]/weigth[vmask].sum())
    Values = []
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(path,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        assert os.path.exists(chirps_file),"the file does not exists {}".format(chirps_file)
        # print(chirps_file)
        data = gdal.Open(chirps_file).ReadAsArray()
        print("indRF",ind)

        pvalues = np.take(data,ind)
        # print(pvalues,weigth/1000000)

        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 :
            print("rf",vmask[0].shape,weigth[vmask].sum() , weigth.sum()/2.0)
            Values.append(defval)
        else:
            wp = np.sum(vweigth*vpvalues)/vweigth.sum()
            Values.append(wp)
        del data
        del pvalues
        del vmask
        del vweigth
        del vpvalues

    return Values

def ExtractPVI(feature,monthtype,defval=-9999):
    path = r"D:\Cornell\EthiopianDrought\AData\PVI10day"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI10day\long_pvi_2003.tif")
    geot = raster.GetGeoTransform()
    XSize, YSize = raster.RasterXSize, raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src, proj)
    invct = osr.CoordinateTransformation(proj, src)
    Centroid = [feature.Centroid().GetX(), feature.Centroid().GetY()]
    centroidY, centroidX = ct.TransformPoint(Centroid[0], Centroid[1])[0:2]
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
    print("are", feature.Area() / 1000000, ccpixelpoly.Area() / 1000000)
    if feature.Area() < 0.2 * ccpixelpoly.Area():
        return None

    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1] - geot[0]) / geot[1] for p in NewPoints]
    Row = [(p[0] - geot[3]) / geot[5] for p in NewPoints]

    minCol, maxCol = np.array(Col).min(), np.array(Col).max()
    minRow, maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol, maxCol = int(minCol), int(np.ceil(maxCol))
    minRow, maxRow = int(minRow), int(np.ceil(maxRow))

    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol, maxCol + 1):
        for irow in range(minRow, maxRow + 1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pc1, pr1 = geot[0] + icol * geot[1], geot[3] + irow * geot[5]
            pc2, pr2 = geot[0] + (icol + 1) * geot[1], geot[3] + irow * geot[5]
            pc3, pr3 = geot[0] + (icol + 1) * geot[1], geot[3] + (irow + 1) * geot[5]
            pc4, pr4 = geot[0] + icol * geot[1], geot[3] + (irow + 1) * geot[5]

            pjc1, pjr1 = invct.TransformPoint(pr1, pc1)[0:2]
            pjc2, pjr2 = invct.TransformPoint(pr2, pc2)[0:2]
            pjc3, pjr3 = invct.TransformPoint(pr3, pc3)[0:2]
            pjc4, pjr4 = invct.TransformPoint(pr4, pc4)[0:2]

            pixel.AddPoint(pjc1, pjr1)
            pixel.AddPoint(pjc2, pjr2)
            pixel.AddPoint(pjc3, pjr3)
            pixel.AddPoint(pjc4, pjr4)
            pixel.AddPoint(pjc1, pjr1)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if feature.Intersect(pixelpoly):
                intersectionpoly = feature.Intersection(pixelpoly)
                if intersectionpoly.Area() >= ccpixelpoly.Area() * 0.2:
                    vCol.append(icol)
                    vRow.append(irow)
                    weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()
    if len(vCol) == 0:
        return None
    print(vCol, vRow)
    ind = np.array(vRow) * XSize + np.array(vCol)
    weigth = np.array(weigth)
    # pvalues = np.take(data,ind)

    # print(pvalues)
    # for id,_ in enumerate(vRow):
    #     print()
    #     print(data[vRow[id]][vCol[id]])
    # weigth = np.array(weigth)
    # vmask = np.where(pvalues != defval)
    # print(weigth[vmask],pvalues[vmask],weigth[vmask]*pvalues[vmask]/weigth[vmask].sum())
    Values = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "{}_pvi_{}.tif".format(monthtype,str(dt.year)))
        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        print("ind", ind)

        pvalues = np.take(data, ind)
        print(pvalues, weigth / 1000000)

        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0:
            print("rf", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
            Values.append(defval)
        else:
            wp = np.sum(vweigth * vpvalues) / vweigth.sum()
            Values.append(wp)
        del data
        del pvalues
        del vmask
        del vweigth
        del vpvalues

    return Values
import pandas as pd

def ExtractCSV(path,croptype,CropPolygon):

    keys = [croptype + str(year) for year in range(2010,2017)]
    areakeys = [croptype[:-3] +"AREA" + str(year) for year in range(2010,2017)]
    data = pd.read_csv(os.path.join(path,croptype+"_AREA.csv"))

    CropList = data[keys].to_numpy()
    AreaList = data[areakeys].to_numpy()
    CropId = data["CropID"].to_numpy().tolist()



    CropListV2 = []
    AreaListV2 = []
    CropIdV2 = []
    RFValues = []
    SPVIValues = []
    LPVIValues = []



    shppath = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shppath)
    layer = dataset.GetLayer()

    # 投影坐标系之间的转换
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src, proj)
    invct = osr.CoordinateTransformation(proj, src)


    for index, ID in enumerate(CropId):
        print(ID)

        id = int(ID) - 1
        feature = layer.GetFeature(id)
        Geom = feature.GetGeometryRef()
        fcount = Geom.GetGeometryCount()
        gname = Geom.GetGeometryName()
        Centroid = [Geom.Centroid().GetX(), Geom.Centroid().GetY()]
        centroidY, centroidX = ct.TransformPoint(Centroid[0], Centroid[1])[0:2]


        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(centroidX, centroidY)
        print("croppolygon",CropPolygon)
        print("point",point)


        if not CropPolygon.Contains(point):
            print("skip")
            continue



        if gname == "MULTIPOLYGON":
            continue
        elif gname == "POLYGON" and fcount == 1:

            RF = ExtractRainFall(Geom)
            print("FF")
            if RF == None:
                continue
            RFValues.append(RF)
            SPVI = ExtractPVI(Geom,"short")

            if SPVI == None:
                break
            SPVIValues.append(SPVI)

            LPVI = ExtractPVI(Geom, "long")
            if LPVI == None:
                break
            LPVIValues.append(LPVI)
            CropListV2.append(CropList[index,:].tolist())
            AreaListV2.append(AreaList[index, :].tolist())
            CropIdV2.append(CropId[index])




    CropListV2 = np.array(CropListV2)
    AreaListV2 = np.array(AreaListV2)
    CropIdV2 = np.array(CropIdV2)
    RFValues = np.array(RFValues)
    SPVIValues = np.array(SPVIValues)
    LPVIValues = np.array(LPVIValues)


    Month =[1,2,3,4,5,6,7,8,9,10,11,12]
    Year =[2010,2011,2012,2013,2014,2015,2016]
    YM = []
    for year in Year:
        for m in Month:
            YM.append(str(year)+str(m).zfill(2))
    DataDict = {}

    for year in range(2010,2017):
        DataDict[croptype+str(year)] = CropListV2[:,year-2010]
        DataDict[croptype[:-3]+"AREA" + str(year)] = AreaListV2[:, year - 2010]

    DataDict["ID"] = CropIdV2

    for idx, ym in enumerate(YM):
        DataDict["RF"+ym] = RFValues[:,idx]


    for idx, ym in enumerate(Year):
        DataDict["ShortPVI"+str(ym)] = SPVIValues[:,idx]
        DataDict["LongPVI" + str(ym)] = LPVIValues[:, idx]


    df = pd.DataFrame(DataDict)
    outpath = r"D:\Cornell\EthiopianDrought\CropCSV\PolyGonCase\{}_VegIndex_Big.csv".format(croptype)
    df.to_csv(outpath, index=False)

croppath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop"



CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH","CropAve"]

daShapefile = r"D:\Cornell\EthiopianDrought\westhihland\westhighland.shp"
# daShapefile = r"D:\Cornell\EthiopianDrought\northhighlandshp\NorthHighland.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()

for ctp in CTp:
    ExtractCSV(croppath,ctp,geo)





