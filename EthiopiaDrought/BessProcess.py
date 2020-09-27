from osgeo import gdal,osr,ogr
from datetime import *
from dateutil import rrule
import os
from epdatapro import write_Img
# import datetime
import glob
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import calendar
def dayofyear(y,m,d):
    """
    :param y:
    :param m:
    :param d:
    :return: the nth day in the year
    """
    days_in_the_year = (date(y, m, d) - date(y, 1, 1)).days + 1
    return str(y)+str(days_in_the_year).zfill(3)
def monthlyComp(directory):
    start = datetime.strptime("-".join(["2000", '01', "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):

        mrange = calendar.monthrange(dt.year,dt.month)[1]
        mstart = dt.date()
        mend = mstart + timedelta(days=mrange-1)
        extendmstart = mstart - timedelta(days=7)
        print(dt)
        vpath = []
        for vdt in (rrule.rrule(rrule.DAILY,interval=1,dtstart=extendmstart,until=mend)):

            if os.path.exists(os.path.join(directory,str(vdt.year)+str(vdt.month).zfill(2)+str(vdt.day).zfill(2)+'.tif')):
                temp =[]
                temp.append(os.path.join(directory,str(vdt.year)+str(vdt.month).zfill(2)+str(vdt.day).zfill(2)+'.tif'))
                if vdt.date() <= mstart:
                    temp.append(8-(mstart-vdt.date()).days)
                elif (vdt+timedelta(days=7)).date() >= mend:
                    temp.append((mend-vdt.date()).days+1)
                else:
                    temp.append(8)
                vpath.append(temp)
        if len(vpath) > 3:
            MosaicImg = np.full(shape=(2400, 2400,len(vpath)), fill_value=255.0)
            Mask = np.full(shape=(2400, 2400,len(vpath)), fill_value=1.0)
            MaskNum = np.full(shape=(2400, 2400, len(vpath)), fill_value=1.0)
            Weigth = np.array([val[1] for val in vpath])
            Weigth = Weigth/Weigth.sum()
            for band in range(len(vpath)):
                Mask[:,:,band] = Mask[:,:,band]*Weigth[band]
            print(Weigth)
            break

            # for band,path in enumerate(vpath):
            #     # 将每个数据插入MosaicImg中
            #     MosaicImg[:,:,band] = gdal.Open(path[0]).ReadAsArray()
            #     # 如果某个像素为nodata，则Mask设置为0，MaskNum也设置为0，MaskNum为了标识有效观测数据个数
            #     nodataMask = np.where(MosaicImg[:,:,band] == 255)
            #     Mask[:,:,band][nodataMask] = 0
            #     MaskNum[:,:,band][nodataMask] = 0
            #
            # for band in range(len(vpath)):
            #     # 加权 ，其中nodata数据位置，为0
            #     MosaicImg[:, :, band] = MosaicImg[:, :, band] * Mask[:, :, band]
            #
            # # 将某一位置所有像素加起来，权重之和加起来，观测个数加起来，
            # MosaicImgSum = MosaicImg.sum(axis=2)
            # MaskWeigthSum = Mask.sum(axis=2)
            # MaskSum = MaskNum.sum(axis=2)
            #
            # # 获取观测数量不够的像素位置
            # noMask = np.where(MaskSum < 3)
            # MosaicImgSum[noMask] = 255
            # MaskWeigthSum[noMask] = 1
            # MosaicImgMean = MosaicImgSum/MaskWeigthSum
            #
            # preference = gdal.Open(r'D:\Cornell\MOD13A3\20000201.tif')
            # proj = osr.SpatialReference()
            # proj.ImportFromWkt(str(preference.GetProjection()))
            # geotrans = preference.GetGeoTransform()
            # path = os.path.join(r'D:\Cornell\BessMonth', str(dt.year) + str(dt.month).zfill(2) + str(1).zfill(2) + '.tif')
            # print(path)
            # write_Img(MosaicImgMean, path, proj, geotrans, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)
            # del MosaicImg
            # del MosaicImgSum
            # del MaskWeigthSum
            # del MaskSum
            # del noMask
            # del MosaicImgMean


# monthlyComp(r'D:\Cornell\Bess')

def mosaic(directory):

    for year in range(2000,2016):
        start = datetime.strptime('{}-01-01'.format(year),'%Y-%m-%d').date()

        stop = datetime.strptime('{}-12-31'.format(year),'%Y-%m-%d').date()

        for dt in (rrule.rrule(rrule.DAILY,interval=8,dtstart=start,until=stop)):
            mdayofyear = dayofyear(dt.year, dt.month, dt.day)
            h21v07 = glob.glob(os.path.join(directory,'BESSET.A{}.h21v07.001.h5'.format(mdayofyear)))
            h21v08 = glob.glob(os.path.join(directory,'BESSET.A{}.h21v08.001.h5'.format(mdayofyear)))
            h22v07 = glob.glob(os.path.join(directory,'BESSET.A{}.h22v07.001.h5'.format(mdayofyear)))
            h22v08 = glob.glob(os.path.join(directory,'BESSET.A{}.h22v08.001.h5'.format(mdayofyear)))

            if (dt.year == 2000 and dt.month == 1) or (dt.year == 2000 and dt.month == 2 and dt.day <26):
                continue
            print(dt.year,dt.month,dt.day)
            assert len(h21v07)*len(h21v08)*len(h22v07)*len(h22v08) >0,"no such file {}".format(dt)

            h21v07et = gdal.Open(h21v07[0]).ReadAsArray().T
            h21v08et = gdal.Open(h21v08[0]).ReadAsArray().T
            h22v07et = gdal.Open(h22v07[0]).ReadAsArray().T
            h22v08et = gdal.Open(h22v08[0]).ReadAsArray().T

            MosaicImg = np.full(shape=(2400,2400),fill_value=255)
            MosaicImg[0:1200,0:1200] = h21v07et
            MosaicImg[1200:2400,0:1200] = h21v08et
            MosaicImg[0:1200, 1200:2400] = h22v07et
            MosaicImg[1200:2400,1200:2400] = h22v08et
        #
            preference = gdal.Open(r'D:\Cornell\MOD13A3\20000201.tif')
            proj = osr.SpatialReference()
            proj.ImportFromWkt(str(preference.GetProjection()))
            geotrans = preference.GetGeoTransform()
            path = os.path.join(directory,str(dt.year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)+'.tif')

            write_Img(MosaicImg, path, proj, geotrans, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)
            del MosaicImg
            del h22v08et
            del h21v08et
            del h22v07et
            del h21v07et



#
# mosaic(r'D:\Cornell\Bess')
# print("I am done")
def BessAggCSV(input):

    minLat,maxLat,minLon,maxLon = [2.5,15,32,49.5]

    latVec = np.arange(maxLat,minLat,-0.25)
    lonVec = np.arange(minLon,maxLon,0.25)
    # print(latVec,len(latVec))
    # print(lonVec,len(lonVec))
    EthImg = np.full(shape=(50,70),fill_value=255.0,dtype=np.float)

    ET = gdal.Open(input)
    ETarr = ET.ReadAsArray()
    BessProjection = ET.GetProjection()
    BessGeost = ET.GetGeoTransform()
    # print(BessGeost)
    gleamPrj = osr.SpatialReference()
    gleamPrj.ImportFromEPSG(4326)
    BessPrj = osr.SpatialReference()
    BessPrj.ImportFromWkt(str(BessProjection))
    geost =[32,0.25,0,15,0,-0.25]
    ct = osr.CoordinateTransformation(gleamPrj,BessPrj)

    ArowList = []
    AcolList = []
    ARList = []
    ACList = []


    for R,lat in enumerate(latVec[1:2]):
        for COL,lon in enumerate(lonVec[1:2]):

            p1lon, p1lat = [lon,lat]
            p2lon, p2lat = [lon+0.25,lat]
            p3lon, p3lat = [lon+0.25,lat-0.25]
            p4lon, p4lat = [lon,lat-0.25]

            p1x, p1y = ct.TransformPoint(p1lat, p1lon)[0:2]
            p2x, p2y = ct.TransformPoint(p2lat, p2lon)[0:2]
            p3x, p3y = ct.TransformPoint(p3lat, p3lon)[0:2]
            p4x, p4y = ct.TransformPoint(p4lat, p4lon)[0:2]

            # print(p1x, p1y,p2x, p2y,p3x, p3y, p4x, p4y)

            p1r, p1c = (p1y - BessGeost[3]) / BessGeost[5], (p1x - BessGeost[0]) / BessGeost[1]
            p2r, p2c = (p2y - BessGeost[3]) / BessGeost[5], (p2x - BessGeost[0]) / BessGeost[1]
            p3r, p3c = (p3y - BessGeost[3]) / BessGeost[5], (p3x - BessGeost[0]) / BessGeost[1]
            p4r, p4c = (p4y - BessGeost[3]) / BessGeost[5], (p4x - BessGeost[0]) / BessGeost[1]

            # print(p1r, p1c, p2r, p2c, p3r, p3c, p4r, p4c)
            minr,maxr = min(p1r,p2r,p3r,p4r) , max(p1r,p2r,p3r,p4r)
            minc,maxc = min(p1c,p2c,p3c,p4c) , max(p1c,p2c,p3c,p4c)

            minr, maxr = int(np.floor(minr)),int(np.ceil(maxr))
            minc, maxc = int(np.floor(minc)),int(np.ceil(maxc))

            plt.imshow(ETarr)
            # plt.plot([p1c,p2c,p3c,p4c,p1c],[p1r,p2r,p3r,p4r,p1r])
            # print(minr,maxr,minc,maxc)
            # Create ring
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(p1r, p1c)
            ring.AddPoint(p4r, p4c)
            ring.AddPoint(p3r, p3c)
            ring.AddPoint(p2r, p2c)
            ring.AddPoint(p1r, p1c)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            plt.plot([p1c,p2c,p3c,p4c,p1c],[p1r,p2r,p3r,p4r,p1r])
            rowList = []
            colList = []
            RList =   []
            CList =   []

            err = 1e-3
            for r in range(minr,maxr+1):
                for col in range(minc,maxc+1):

                    r1, c1 = r, col
                    r2, c2 = r, col + 1 -err
                    r3, c3 = r + 1 -err, col + 1 -err
                    r4, c4 = r +1 - err, col
                    # Create ring
                    mring = ogr.Geometry(ogr.wkbLinearRing)
                    mring.AddPoint(r1, c1)
                    mring.AddPoint(r4, c4)
                    mring.AddPoint(r3, c3)
                    mring.AddPoint(r2, c2)
                    mring.AddPoint(r1, c1)
                    mpoly = ogr.Geometry(ogr.wkbPolygon)
                    mpoly.AddGeometry(mring)
                    if mpoly.Intersect(poly):
                        rowList.append(r1)
                        colList.append(c1)
            plt.scatter(colList,rowList)
            hh = ETarr[(np.array(rowList),np.array(colList))]
            print(list(zip(rowList,colList,hh)))
            print(hh[hh < 255].mean())
            plt.show()
            if len(rowList) > 0:
                RList.append(R)
                CList.append(COL)
                ArowList.append(rowList)
                AcolList.append(colList)
                ARList.append(RList)
                ACList.append(CList)
    # return ArowList,AcolList,ARList,ACList
BessAggCSV(r"D:\Cornell\BessMonth\20000301.tif")




def BessAgg(inputDir,outDir):

    minLat,maxLat,minLon,maxLon = [2.5,15,32,49.5]
    ArowList, AcolList, ARList, ACList = BessAggCSV(r"D:\Cornell\BessMonth\20000301.tif")
    print("start")
    gleamPrj = osr.SpatialReference()
    gleamPrj.ImportFromEPSG(4326)
    geost = [32, 0.25, 0, 15, 0, -0.25]

    bessfiles = glob.glob(os.path.join(inputDir,"*.tif"))

    for input in bessfiles:
        print(input)
        EthImg = np.full(shape=(50,70),fill_value=-999.0,dtype=np.float)
        ET = gdal.Open(input)
        ETarr = ET.ReadAsArray()
        for id , rowlist in enumerate(ArowList):

            mask = (np.array(rowlist), np.array(AcolList[id]))
            pvalues = ETarr[mask]
            pvalues=pvalues[pvalues <255]
            # print("p",pvalues.shape)
            if pvalues.shape[0] >400:
                # print("real",[ARList[id]],[ACList[id]])
                EthImg[ARList[id][0]][ACList[id][0]] = pvalues.mean()
        # print(EthImg)
        output = os.path.join(outDir,os.path.basename(input))

        write_Img(EthImg, output, gleamPrj, geost, 70, 50, im_bands=1, dtype=gdal.GDT_Float32)
        print("first{}".format(input))




# BessAgg(r'D:\Cornell\BessMonth',r'D:\Cornell\BessAgg0.25V2')








