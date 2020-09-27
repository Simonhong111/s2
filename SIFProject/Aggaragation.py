from osgeo import  gdal
from glob import glob
import os
from datetime import datetime
import numpy as np


class Sen2Aggragation(object):

    def __init__(self,year,directory):

        self.year = year
        self.directory = directory

    def findValidMonth(self):
        """
        返回给定年份中存在影像数据的月份
        :return:
        """

        #获取哪一年的数据
        filenames = glob(os.path.join(self.directory,self.year+"*.tif"))
        #影像名称
        filebasenames = [os.path.basename(f) for f in filenames]
        #获取哪一年那一月的数据


        # yymmdd_gridID_tileID.tif
        # 20170821_1_T48SWC.tif
        yymmdd = [file.split("_")[0] for file in filebasenames]

        yymmdd = sorted(list(set(yymmdd)))

        print(yymmdd)

        valmonth = [datetime.strptime(myymmdd,"%Y%m%d").month for myymmdd in yymmdd]

        return valmonth

    def aGgaragation(self):

        mMonth = self.findValidMonth()

        for month in mMonth:
            filebymonth = glob(os.path.join(self.directory,self.year+str(month)))



        pass

    def testAgg(self,path):
        pass

    def calMean(self,path):

        # 计算每个波段有效值的均值（此处还未完善）
        raster = gdal.Open(path)
        mask = None
        BandSumArr =[]
        for band in range(raster.RasterCount-1):
            data = raster.GetRasterBand(band+1).ReadAsArray().astype(np.float)
            # data = data * mask
            Sum = np.sum(data)
            BandSumArr.append(Sum)
            data = None


        return BandSumArr


































