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


def ExtractEVI(feature,defval=-3000):
    path = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
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
    # print(vCol, vRow)
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
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        evi_file = os.path.join(path,
                                str(dt.year) + "." + str(dt.month).zfill(2) + ".01" + ".tif")
        assert os.path.exists(evi_file), "the file does not exists {}".format(evi_file)
        data = gdal.Open(evi_file).ReadAsArray()
        print("indEVI", ind)

        pvalues = np.take(data, ind)
        # print("value",pvalues, weigth / 1000000)

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
    # print("final",Values)
    return Values


def ExtractNEWSIF(feature,defval=-9999):
    path = r"D:\Cornell\NewSIF005Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_200208.nc.tif")
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
    # print(vCol, vRow)
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
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        print("indNEWSIF", ind)

        pvalues = np.take(data, ind)
        # print(pvalues, weigth / 1000000)

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



def ExtractGOSIF(feature,defval=32766):
    path = r"D:\Cornell\GOSIFV002Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
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
    # print("are", feature.Area() / 1000000, ccpixelpoly.Area() / 1000000)
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
    # print(vCol, vRow)
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
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "GOSIF_{}.M{}.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        # print("ind", ind)

        pvalues = np.take(data, ind)
        # print(pvalues, weigth / 1000000)

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


def ExtractET(feature,defval=-999):
    path = r"D:\Cornell\EthiopianDrought\AData\ETClip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ETClip\ET.v3.3a198001.tif")
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

    Values = []
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "ET.v3.3a{}{}.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        print("ind***", ind)

        pvalues = np.take(data, ind)
        print("value",pvalues, weigth / 1000000)

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

def ExtractSM(feature,defval=-2):
    path = r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5\ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-199101.tif")
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
    if feature.Area() < 0.4 * ccpixelpoly.Area():
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
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-{}{}.tif".format(str(dt.year),
                                                                                                 str(dt.month).zfill(
                                                                                                     2)))

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
def ExtractPVI(feature,defval=-9999):
    path = r"D:\Cornell\EthiopianDrought\AData\PVI"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI\pvi_2010.tif")
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
        newsif_file = os.path.join(path, "long_pvi_{}.tif".format(str(dt.year)))
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

def ExtractCSV(path):
    CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]
    CropType = []
    for mctp in CTp:
        for year in range(2010, 2017):
            CropType.append(mctp + str(year))
            CropType.append(mctp[:-3] + "AREA" + str(year))

    CropList = []
    CropId = []

    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            temp = []
            for mcp in CropType:
                temp.append(row[mcp])
            CropId.append(row["ID"])
            CropList.append(temp)

    CropListV2 = []
    CropIdV2 = []
    RFValues = []
    EVIValues = []
    NEWSIFValues = []
    GOSIFValues = []
    ETValues = []
    SMValues = []
    PVIValues = []

    shppath = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shppath)
    layer = dataset.GetLayer()

    for index, ID in enumerate(CropId):
        print(ID)

        id = int(ID) - 1
        feature = layer.GetFeature(id)
        Geom = feature.GetGeometryRef()
        fcount = Geom.GetGeometryCount()
        gname = Geom.GetGeometryName()
        # print("geometrytype", gname)

        if gname == "MULTIPOLYGON":
            continue
        elif gname == "POLYGON" and fcount == 1:

            # RF = ExtractRainFall(Geom)
            # print("FF")
            # if RF == None:
            #     continue
            # RFValues.append(RF)
            EVI = ExtractEVI(Geom)
            if EVI == None:
                continue
                # break

            EVIValues.append(EVI)
            # NSIF = ExtractNEWSIF(Geom)
            # if NSIF == None:
            #     break
            #
            # NEWSIFValues.append(NSIF)
            # GSIF = ExtractGOSIF(Geom)
            # if GSIF == None:
            #     break
            # GOSIFValues.append(GSIF)
            # ET = ExtractET(Geom)
            # if ET == None:
            #     break
            # ETValues.append(ET)
            # SM = ExtractSM(Geom)
            # if SM == None:
            #     break
            # SMValues.append(SM)
            # PVI = ExtractPVI(Geom)
            # if PVI == None:
            #     break
            # PVIValues.append(PVI)
            # CropListV2.append(CropList[index])
            # CropIdV2.append(CropId[index])


    # VValues = []
    # Weigth = []
    # for ccp in CropListV2:
    #
    #     VValues.append([ccp[i*2] for i in range(35)])
    #     Weigth.append([ccp[i*2+1] for i in range(35)])
    #
    #
    #
    # VValues = np.array(VValues).astype(np.float)
    # Weigth = np.array(Weigth).astype(np.float)
    # mCrop2010 = VValues[:,[7*i for i in range(5)]]
    # mWeigth2010 = Weigth[:,[7*i for i in range(5)]]
    # mCrop2010 = np.sum(mCrop2010*mWeigth2010,axis=1)/np.sum(mWeigth2010,axis=1)
    #
    # mCrop2011 = VValues[:,[7*i+1 for i in range(5)]]
    # mWeigth2011 = Weigth[:,[7*i+1 for i in range(5)]]
    # mCrop2011 = np.sum(mCrop2011 * mWeigth2011, axis=1) / np.sum(mWeigth2011, axis=1)
    #
    # mCrop2012 = VValues[:,[7*i+2 for i in range(5)]]
    # mWeigth2012 = Weigth[:,[7*i+2 for i in range(5)]]
    # mCrop2012 = np.sum(mCrop2012 * mWeigth2012, axis=1) / np.sum(mWeigth2012, axis=1)
    #
    # mCrop2013 = VValues[:,[7*i+3 for i in range(5)]]
    # mWeigth2013 = Weigth[:,[7*i+3 for i in range(5)]]
    # mCrop2013 = np.sum(mCrop2013 * mWeigth2013, axis=1) / np.sum(mWeigth2013, axis=1)
    #
    # mCrop2014 = VValues[:,[7*i+4 for i in range(5)]]
    # mWeigth2014 = Weigth[:,[7*i+4 for i in range(5)]]
    # mCrop2014 = np.sum(mCrop2014 * mWeigth2010, axis=1) / np.sum(mWeigth2014, axis=1)
    #
    # mCrop2015 = VValues[:,[7*i+5 for i in range(5)]]
    # mWeigth2015 = Weigth[:,[7*i+5 for i in range(5)]]
    # mCrop2015 = np.sum(mCrop2015 * mWeigth2015, axis=1) / np.sum(mWeigth2015, axis=1)
    #
    # mCrop2016 = VValues[:,[7*i+6 for i in range(5)]]
    # mWeigth2016 = Weigth[:,[7*i+6 for i in range(5)]]
    # mCrop2016 = np.sum(mCrop2016 * mWeigth2016, axis=1) / np.sum(mWeigth2016, axis=1)
    #
    #
    #
    #
    #
    # CropListV2 = np.array(CropListV2)
    # CropIdV2 = np.array(CropIdV2)
    # RFValues = np.array(RFValues)
    # EVIValues = np.array(EVIValues)
    # NEWSIFValues = np.array(NEWSIFValues)
    # GOSIFValues = np.array(GOSIFValues)
    # # ETValues = np.array(ETValues)
    # # SMValues = np.array(SMValues)
    # PVIValues = np.array(PVIValues)
    #
    #


    # Month =[1,2,3,4,5,6,7,8,9,10,11,12]
    # Year =[2010,2011,2012,2013,2014,2015,2016]
    # YM = []
    # for year in Year:
    #     for m in Month:
    #         YM.append(str(year)+str(m).zfill(2))
    # DataDict = {"mCrop2010":mCrop2010,"mCrop2011":mCrop2011,"mCrop2012":mCrop2012,"mCrop2013":mCrop2013,
    #             "mCrop2014":mCrop2014,
    #             "mCrop2015":mCrop2015,"mCrop2016":mCrop2016}
    #
    #
    #
    # for idx,croptype2 in enumerate(CropType):
    #     DataDict[croptype2] = CropListV2[:,idx]
    # DataDict["ID"] = CropIdV2
    #
    # for idx, ym in enumerate(YM):
    #     DataDict["RF"+ym] = RFValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["EVI"+ym] = EVIValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["NEWSIF"+ym] = NEWSIFValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["GOSIF"+ym] = GOSIFValues[:,idx]
    # # for idx, ym in enumerate(YM):
    # #     print("etva",ETValues)
    # #     DataDict["ET"+ym] = ETValues[:,idx]
    # # for idx, ym in enumerate(YM):
    # #     DataDict["SM"+ym] = SMValues[:,idx]
    #
    # for idx, ym in enumerate(Year):
    #
    #     DataDict["PVI"+str(ym)] = PVIValues[:,idx]
    # print(DataDict["PVI2010"])
    # df = pd.DataFrame(DataDict)
    # outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\CROP_Index_From2010-2016_V2long.csv"
    # df.to_csv(outpath, index=False)

croppath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Crop_2010-2016.csv"


ExtractCSV(croppath)



