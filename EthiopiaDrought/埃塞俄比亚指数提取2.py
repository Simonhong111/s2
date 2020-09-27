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
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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

        data = gdal.Open(chirps_file).ReadAsArray()
        pvalues = np.take(data,ind)

        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
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

    # for i in range(len(vRow)):
    #     plt.plot(vCol[i],vRow[i])
    # colrow = [156.0,154.35668307791,156.003698973185,154.358221310937,156.039739094546,154.406611900221,156.058527268024,
    #           154.454935206114,156.122358118497,154.455184435339,156.143490106408,154.45245409785,156.14824896573,154.45183922251,
    #           156.160602157942,154.381235947108,156.160690048567,154.358834552837,156.195498115771,154.281426516919,156.209496080576,
    #           154.231508410439,156.221714814625,154.195369013248,156.218474624376,154.141936895723,156.201371692443,154.10395919113,
    #           156.175615393369,154.072840356093,156.170690514732,154.0090623215,156.16891293433,154.0,156,154,156.0,154.35668307791]
    # colrow = np.array(colrow)
    # plt.plot(colrow[::2],colrow[1::2],"k-")
    # plt.show()

def ExtractEVI(feature,defval=-3000):
    path = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
        evi_file = os.path.join(path,
                                str(dt.year) + "." + str(dt.month).zfill(2) + ".01" + ".tif")
        assert os.path.exists(evi_file), "the file does not exists {}".format(evi_file)
        data = gdal.Open(evi_file).ReadAsArray()
        pvalues = np.take(data,ind)

        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("evi", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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
def ExtractNEWSIF(feature,defval=-9999):
    path = r"D:\Cornell\NewSIF005Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_200208.nc.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
        newsif_file = os.path.join(path, "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        pvalues = np.take(data,ind)
        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]

        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("newsif", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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
def ExtractGOSIF(feature,defval=32766):
    path = r"D:\Cornell\GOSIFV002Clip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
        newsif_file = os.path.join(path, "GOSIF_{}.M{}.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        pvalues = np.take(data,ind)
        vmask = np.where(pvalues < defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]

        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("gosif", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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
def ExtractET(feature,defval=-999):
    path = r"D:\Cornell\EthiopianDrought\AData\ETClip"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ETClip\ET.v3.3a198001.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
        newsif_file = os.path.join(path, "ET.v3.3a{}{}.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        pvalues = np.take(data,ind)
        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("et", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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
def ExtractSM(feature,defval=-2):
    path = r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5\ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-199101.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
        newsif_file = os.path.join(path, "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-{}{}.tif".format(str(dt.year),
                                                                                                 str(dt.month).zfill(
                                                                                                     2)))

        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        pvalues = np.take(data,ind)
        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("et", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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
def ExtractPVI(feature,defval=-9999):
    path = r"D:\Cornell\EthiopianDrought\AData\PVI"
    start = datetime.strptime("-".join(["2010", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2016-12-31", "%Y-%m-%d").date()
    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI\pvi_2010.tif")
    geot = raster.GetGeoTransform()
    XSize,YSize = raster.RasterXSize,raster.RasterYSize
    proj = osr.SpatialReference()
    proj.ImportFromWkt(str(raster.GetProjection()))
    src = osr.SpatialReference()
    src.ImportFromEPSG(32637)
    ct = osr.CoordinateTransformation(src,proj)
    Points = Poy2Points(feature)
    NewPoints = ct.TransformPoints(Points)
    Col = [(p[1]-geot[0])/geot[1] for p in NewPoints]
    Row = [(p[0]-geot[3])/geot[5] for p in NewPoints]

    minCol,maxCol = np.array(Col).min(), np.array(Col).max()
    minRow,maxRow = np.array(Row).min(), np.array(Row).max()
    # print(minCol, maxCol, minRow, maxRow)
    minCol,maxCol = int(minCol),int(np.ceil(maxCol))
    minRow,maxRow = int(minRow),int(np.ceil(maxRow))
    # print(minCol,maxCol,minRow,maxRow)
    # plt.imshow(raster.ReadAsArray())
    # plt.plot(Col,Row)
    # plt.show()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id,_ in enumerate(Row):
        ring.AddPoint(Col[id],Row[id])
    ring.AddPoint(Col[0],Row[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    # area = poly.Area()
    # print(area)
    vCol = []
    vRow = []
    weigth = []
    for icol in range(minCol,maxCol+1):
        for irow in range(minRow,maxRow+1):
            pixel = ogr.Geometry(ogr.wkbLinearRing)
            pixel.AddPoint(icol,irow)
            pixel.AddPoint(icol+1, irow)
            pixel.AddPoint(icol+1, irow+1)
            pixel.AddPoint(icol , irow+1)
            pixel.AddPoint(icol , irow)
            pixelpoly = ogr.Geometry(ogr.wkbPolygon)
            pixelpoly.AddGeometry(pixel)
            if poly.Intersect(pixelpoly):
                # vCol.append([icol,icol+1,icol+1,icol,icol])
                # vRow.append([irow,irow,irow+1,irow+1,irow])
                vCol.append(icol)
                vRow.append(irow)
                intersectionpoly = poly.Intersection(pixelpoly)
                weigth.append(intersectionpoly.Area())
    # data = raster.ReadAsArray()

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
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        newsif_file = os.path.join(path, "pvi_{}.tif".format(str(dt.year)))
        assert os.path.exists(newsif_file), "the file does not exists {}".format(newsif_file)
        data = gdal.Open(newsif_file).ReadAsArray()
        pvalues = np.take(data,ind)
        vmask = np.where(pvalues != defval)
        vweigth = weigth[vmask]
        vpvalues = pvalues[vmask]
        print(dt, newsif_file)
        print(vCol, vRow, weigth, pvalues, np.sum(pvalues * weigth) / weigth.sum())
        # print("shape",vmask[0].shape[0])
        if vmask[0].shape[0] == 0 or weigth[vmask].sum() < weigth.sum()/2.0:
            print("pvi", vmask[0].shape, weigth[vmask].sum(), weigth.sum() / 2.0)
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

# aggpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop"
# CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
# crtype = CTp[0]
# aggtypepath = os.path.join(aggpath,crtype+".csv")
# CropList = []
# CropId = []
# CropType = os.path.basename(aggtypepath)[:-4]
# with open(aggtypepath, newline='') as csvfile:
#     reader = csv.DictReader(csvfile)
#     for row in reader:
#         temp = []
#         temp.append(row[CropType+"2010"])
#         temp.append(row[CropType+"2011"])
#         temp.append(row[CropType+"2012"])
#         temp.append(row[CropType+"2013"])
#         temp.append(row[CropType+"2014"])
#         temp.append(row[CropType+"2015"])
#         temp.append(row[CropType+"2016"])
#         CropId.append(row["CropID"])
#         CropList.append(temp)
#
# shppath = r"D:\Cornell\EthiopianDrought\CropCSV\sub_kebele_shapefiles\Export_Output.shp"
# driver = ogr.GetDriverByName("ESRI Shapefile")
# dataset = driver.Open(shppath)
# layer = dataset.GetLayer()
#
# CropListV2 = []
# CropIdV2 = []
# for index,ID in enumerate(CropId):
#
#     id = int(ID) - 1
#     feature = layer.GetFeature(id)
#     Geom = feature.GetGeometryRef()
#     fcount = Geom.GetGeometryCount()
#     gname = Geom.GetGeometryName()
#     print("geometrytype", gname)
#
#     if gname == "MULTIPOLYGON":
#         continue
#     elif gname == "POLYGON" and fcount ==1:
#         ExtractRainFall(Geom)

def ExtractCSV(path):

    CropList = []
    CropId = []
    CropType = os.path.basename(path)[:-4]
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
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
            CropListV2.append(CropList[index])
            CropIdV2.append(CropId[index])
            # RFValues.append(ExtractRainFall(Geom))
            # EVIValues.append(ExtractEVI(Geom))
            # NEWSIFValues.append(ExtractNEWSIF(Geom))
            # GOSIFValues.append(ExtractGOSIF(Geom))
    #         ETValues.append(ExtractET(Geom))
    #         SMValues.append(ExtractSM(Geom))
            PVIValues.append(ExtractPVI(Geom))
    #
    #
    # CropListV2 = np.array(CropListV2)
    # CropIdV2 = np.array(CropIdV2)
    # RFValues = np.array(RFValues)
    # EVIValues = np.array(EVIValues)
    # NEWSIFValues = np.array(NEWSIFValues)
    # GOSIFValues = np.array(GOSIFValues)
    # ETValues = np.array(ETValues)
    # SMValues = np.array(SMValues)
    # PVIValues = np.array(PVIValues)
    # Month =[1,2,3,4,5,6,7,8,9,10,11,12]
    # Year =[2010,2011,2012,2013,2014,2015,2016]
    # YM = []
    # for year in Year:
    #     for m in Month:
    #         YM.append(str(year)+str(m).zfill(2))
    # DataDict = {CropType + '2010': CropListV2[:, 0], CropType + '2011': CropListV2[:, 1], CropType + '2012': CropListV2[:, 2],
    #       CropType + '2013': CropListV2[:, 3], CropType + '2014': CropListV2[:, 4], CropType + '2015': CropListV2[:, 5],
    #       CropType + '2016': CropListV2[:, 6], "CropID": CropIdV2}
    # for idx, ym in enumerate(YM):
    #     DataDict["RF"+ym] = RFValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["EVI"+ym] = EVIValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["NEWSIF"+ym] = NEWSIFValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["GOSIF"+ym] = GOSIFValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["ET"+ym] = ETValues[:,idx]
    # for idx, ym in enumerate(YM):
    #     DataDict["SM"+ym] = SMValues[:,idx]
    # for idx, ym in enumerate(Year):
    #     DataDict["PVI"+str(ym)] = PVIValues[:,idx]
    #
    # df = pd.DataFrame(DataDict)
    # outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\{}_From2010-2016_V2.csv".format(CropType)
    # df.to_csv(outpath, index=False)
CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
croppath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop"

for ctp in CTp[0:1]:
    mpath = os.path.join(croppath,ctp+'.csv')
    ExtractCSV(mpath)



