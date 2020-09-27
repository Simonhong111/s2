#!/usr/bin/env python
#
# Copyright 2019 
# This program is free software: you can redistribute it and/or modify
# it
# Author: Simon Hong & Yineasttoei Zhang
# Datetime May 1st 2019


from osgeo import gdal, osr
import glob, os, sys
import numpy as np
from geo2mapxy import *


class SenL2AReader(object):
    def __init__(self,Sen2Directory):
        """
        :param path: 指向L2A 文件的路径
        :param res: L2A 数据分辨率
        :param qi:

        """
        self.Directory = Sen2Directory
        self.RasterRes = 20
        self.QiRes = 20
        self.QiDataPath = self.getQiDataPath()
        self.ImgDataPath = self.getImgDataPath()
        self.RasterData = None
        self.QiData = None
        self.RasterXsize = 0
        self.RasterYsize = 0
        self.GeoTransForm = None
        self.ProjectionRef = None

    def setRasterResolution(self,RsResolution):
        """
        设置栅格数据分辨率
        :param RsResolution:
        :return:
        """


        self.RasterRes = RsResolution

    def setQiResolution(self,QiResolution):

        """
        设置质量控制图像的分辨率
        :param QiResolution:
        :return:
        """
        self.QiRes = QiResolution

    def setClassParameter(self):
        """
        设置公用参数，也可以单独设置
        :return:
        """
        if self.RasterData is None:
            raise Exception("The RasterData is None, please call setRasterData")
        else:
            self.setRasterXYSize()
            self.setGeoTransForm()
            self.setProjectionRef()



    def getImgDataPath(self):
        """
        获取大气校正后 IMG_DATA 数据路径
        :return:
        """

        sen2datapath = os.path.join(glob.glob(os.path.join(self.Directory, "GRANULE", "L2*"))[0],"IMG_DATA")
        return sen2datapath

    def getQiDataPath(self):
        """
        获取 QI_DATA 质量控制数据路径
        :return:
        """
        sen2qidatapath = os.path.join(glob.glob(os.path.join(self.Directory, "GRANULE", "L2*"))[0], "QI_DATA")
        return sen2qidatapath

    def setRasterData(self,filename):
        """
        设置公用栅格数据，方便后续读取数据元数据
        :param filename: 栅格数据名称（B02，B03，...）
        :return:
        """

        self.RasterData = self.getImgData(filename)

    def getImgData(self,filename):
        """
        获取IMG_DATA 数据
        :param filename:
        :return:
        """
        # 定位到 raster 数据路径
        respath = os.path.join(self.ImgDataPath,"R"+str(self.RasterRes)+"m")

        # 定位到具体波段，如 B02 波段
        rasterpath = glob.glob(os.path.join(respath,'*'+filename+"*"))

        # 判断所给定的波段名称是否存在
        if len(rasterpath) != 1:
            raise Exception("The band name is wrong, please give correct band name such as B02,B03,...,B8A")

        try:
            raster = gdal.Open(rasterpath[0])
        except Exception as e:
            print("The {} cannot be read in successfully".format(filename))
            print(e)

        return raster


    def setQiData(self,filename):
        """
        设置质量控制数据
        :param filename:质量控制文件名称
        :return:
        """

        self.QiData = self.getQiData(filename)

    def setQiDataFromImgData(self,filename):

        self.QiData = self.getImgData(filename)

    def getQiData(self,filename):
        """
        获取QI_DATA 质量控制数据
        :param filename:
        :return:
        """
        # 定位到质量文件
        qidatapath = glob.glob(os.path.join(self.QiDataPath, '*' + filename + "_" + str(self.QiRes) + "*"))

        # 判断质量控制数据是否存在
        if len(qidatapath) != 1:
            raise Exception("The qi data name is wrong, please give correct band name such as SLC,CLDPRB")

        try:
            qidata = gdal.Open(qidatapath[0])
        except Exception as e:
            print(e)

        return qidata

    def setGeoTransForm(self):

        try:
            geotransform = self.RasterData.GetGeoTransform()
        except Exception as e:
            print("Maybe The RasterData is Empty, please setRasterData Firstly")
            raise e

        if geotransform is None or geotransform is []:
            raise Exception("The geotransform is empty")

        self.GeoTransForm = geotransform

    def setProjectionRef(self):
        """
        设置公用的投影坐标信息
        :return:
        """

        try:
            projectionref = self.RasterData.GetProjectionRef()

        except Exception as e:
            print("Maybe The RasterData is Empty, please setRasterData Firstly")
            raise e
        if projectionref is None or projectionref is []:
            raise Exception("The geotransform is empty")

        self.ProjectionRef = projectionref


    # def getRoi(self,outfile,scope):
    #     """
    #     :param outfile:
    #     :param scope:
    #     :return:
    #     """
    #
    #     pszSrcFile = self.RasterData
    #
    #     if pszSrcFile is None:
    #         raise Exception("The source file is empty, please check the self.RasterData ")
    #     else:
    #         pszDstFile = gdal.Translate(outfile, pszSrcFile, projWin=scope)
    #         pszDstFile = None
    #

    def setRasterXYSize(self):
        """
        设置栅格数据的尺度
        :return:
        """
        if self.RasterData is not None:
            self.RasterXsize = self.RasterData.RasterXSize
            self.RasterYsize = self.RasterData.RasterYSize
        else:
            raise Exception("The RasterData is Empty, please call setRasterData firstly")

    def getRasterXYSize(self):
        """
        获得栅格数据的尺度
        :return:
        """
        if self.RasterData is None:
            raise Exception("The RasterData is empty, please call setRasterData firstly")

        return self.RasterXsize,self.RasterYsize


    def isInsectionwithRaster(self, referbox, targetbox):
        """
        判断 targetbox 与 referbox 是否有交集
        :param referbox:
        :param targetbox:
        :return:
        """
        refxcenter = (referbox[0]+referbox[2])/2.0
        referycenter = (referbox[1]+referbox[3])/2.0

        tarxcenter = (targetbox[0]+targetbox[2])/2.0
        tarycenter = (targetbox[1]+targetbox[3])/2.0

        refwidth = np.abs(referbox[0]-refxcenter)
        refheigth = np.abs(referbox[1]-referycenter)

        tarwidth = np.abs(targetbox[0]-tarxcenter)
        tarheigth = np.abs(targetbox[1]-tarycenter)

        centerxdis = np.abs(refxcenter-tarxcenter)
        centerydis = np.abs(referycenter-tarycenter)

        if (refwidth+tarwidth) > centerxdis and (refheigth+tarheigth)> centerydis:
            return True
        else:
            return False


    def isClip(self,bound):
        """
        判断 bound 范围是与raster有重叠，是否可以用bound 裁剪raster
        :param bound: (ulx,uly,lrx,lry)左上角右下角经纬度
        :return: 是否
        """

        uppergeo = lonlat2geo(self.RasterData,bound[0],bound[1])
        bottomgeo = lonlat2geo(self.RasterData,bound[2],bound[3])

        upxy = geo2imagexy(self.RasterData,uppergeo[0],uppergeo[1])
        btxy = geo2imagexy(self.RasterData,bottomgeo[0],bottomgeo[1])

        if self.RasterXsize * self.RasterYsize == 0:
            raise Exception("The rasterX/Ysize is 0, please call setRasterXYSize or setClassParameter ")

        return self.isInsectionwithRaster((0,0,self.RasterXsize,self.RasterYsize),(upxy[0],upxy[1],btxy[0],btxy[1]))

    def maskByQiData(self,criterion,rasterName):
        """
        :param criterion: 质量控制标准
        :param rasterName: 需要掩膜的波段矢量，例如 B02
        :return: 经过质量控制的 Array 数据
        """

        try:

            QiData = self.QiData.ReadAsArray().astype(np.float)
            raster = self.getImgData(rasterName).ReadAsArray().astype(np.float)

        except Exception as e:
            print("RasterData/QiData can not be read in successfully")
            raise e

        if len(criterion) <1:
            print("The criterion is empty, please give the criterion")
            sys.exit()

        mask = (QiData != criterion[0])
        # print("0",mask[1621][3181])
        # print("mask0", mask[0][0])
        # print("mask1", (QiData != criterion[2])[0][0])
        # print( "wht",(QiData != criterion[2])[0][0]^mask[0][0])


        for criter in criterion[1:]:


            mask &= (QiData != criter)




        print("maskend",mask[1621][3181])
        return raster * mask

    def composite(self,outfile, bandlist, criterion, type ="uint16"):
        """
        description: 该函数用于合成多波段数据
        :param outfile: 合成数据输出文件
        :param bandlist: 选择用于合成的波段
        :param criterion: 质量控制掩膜
        :param type: 输出数据类型
        :return:
        """
        bandnum = len(bandlist)

        if self.RasterXsize * self.RasterYsize == 0:

            print("rasterXSize or rasterYSize is 0, please call setRasterXYSize")
            sys.exit()
        else:
            width = self.RasterXsize
            height = self.RasterYsize

        if type in ["int8", "uint8"]:
            datatype = gdal.GDT_Byte
        elif type in ["int16", "uint16"]:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32

        if self.GeoTransForm is None:
            print("GeoTransForm is None, please setGeoTransForm firstly")
            sys.exit()
        else:
            mgeotransform = self.GeoTransForm

        # 影像的 GeoTransform
        # GeoTransForm[0] 影像左上角像元的左上角定点的地理坐标（单位m），东西方向
        # GeoTransForm[3] 影像左上角像元的左上角定点的地理坐标（单位m），南北方向
        # GeoTransForm[1] 影像像元东西方向的宽度，即东西方向分辨率
        # GeoTransForm[5] 影像像元南北方向的宽度，即南北方向分辨率
        # GeoTransForm[2] 影像旋转相关参数，通常为0
        # GeoTransForm[4] 影像旋转相关参数，通常为0

        originX = mgeotransform[0]
        originY = mgeotransform[3]
        pixelWidth = mgeotransform[1]
        pixelHeight = mgeotransform[5]
        # rotate1 = mgeotransform[2]
        # rotate2 = mgeotransform[4]


        driver = gdal.GetDriverByName("GTiff")
        print("driver", driver)

        dataset = driver.Create(outfile, width, height, bandnum, datatype)
        dataset.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

        if dataset is None:
            print("driver is None, please chech it")
            sys.exit()


        for idx,band in enumerate(bandlist):
            print(bandlist[idx])
            mRasterbyMask = self.maskByQiData(criterion,bandlist[idx])

            dataset.GetRasterBand(idx + 1). \
                    WriteArray(mRasterbyMask)
            mRasterbyMask = None

        outRasterSRS = osr.SpatialReference()
        if self.ProjectionRef is None:
            print("projectionref is None, please call setProjectionRed firstly")
            sys.exit()
        outRasterSRS.ImportFromWkt(self.ProjectionRef)
        dataset.SetProjection(outRasterSRS.ExportToWkt())
        print("save ok")




















