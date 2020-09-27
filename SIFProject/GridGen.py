"""
配置湖北省所对应的格网信息（CMG）
# CMG 0.05 x 0.05 degree
#
"""
import os
import numpy as np
from osgeo import gdal,ogr,osr
from GetPointFromShp import *
class GridDefination(object):

    def __init__(self):

        # 最大/最小的经纬度，确定CMG格网划分
        self.MAX_LONGITUDE = 180
        self.MIN_LONGITUDE = -180
        self.MAX_LATITUDE = 90
        self.MIN_LATITUDE = -90

        # 自定义经纬度，确定研究范围在CMG中的格网
        self.USER_MAX_LONGITUDE = 180
        self.USER_MIN_LONGIUDE = -180
        self.USER_MAX_LATITUDE = 90
        self.USER_MIN_LATITUDE = -90

        # 经纬度间隔
        self.LON_INTERVAL =0.05
        self.LAT_INTERVAL =0.05
        self.CMGGRID = None

    def setUserExtent(self,extent):
        """

        :param extent:(ulx,lrx,lry,uly)
        :return:
        """
        self.USER_MAX_LONGITUDE = extent[1]
        self.USER_MIN_LONGIUDE = extent[0]
        self.USER_MAX_LATITUDE = extent[3]
        self.USER_MIN_LATITUDE = extent[2]

    def user_grid(self):
        """
        自定义格网范围
        计算用户自定义范围，对应于CMG格网的范围
        :return: 返回格网配置

        """
        # 网格ID位置
        user_min_Lon_id = int(np.floor((self.USER_MIN_LONGIUDE - self.MIN_LONGITUDE)/self.LON_INTERVAL))
        user_max_Lon_id = int(np.ceil((self.USER_MAX_LONGITUDE - self.MIN_LONGITUDE)/self.LON_INTERVAL))
        user_min_Lat_id = int(np.floor((self.USER_MIN_LATITUDE - self.MIN_LATITUDE) / self.LAT_INTERVAL))
        user_max_Lat_id = int(np.ceil((self.USER_MAX_LATITUDE - self.MIN_LATITUDE) / self.LAT_INTERVAL))



        for lat in range(user_min_Lat_id, user_max_Lat_id):
            for lon in range(user_min_Lon_id, user_max_Lon_id):
                #(最小经度，最小维度，最大经度，最大纬度，小格网id)
                ulx = lon * self.LON_INTERVAL + self.MIN_LONGITUDE
                lrx = ulx + self.LON_INTERVAL
                lry = lat * self.LAT_INTERVAL + self.MIN_LATITUDE
                uly = lry + self.LON_INTERVAL
                # print("ulx,lrx,uly,lry",ulx,lrx,uly,lry)
                extent = yield (ulx,uly,lrx,lry)   # 生成一个方形的里面有好多格网的矩阵 从最后一行开始编号

    def GridToalNum(self):

        """
        定义给定大小的CMG格网
        :return:格网配置，包括范围，经纬度，编号等
        """
        lat_grid_num = (self.MAX_LATITUDE-self.MIN_LATITUDE)/self.LAT_INTERVAL
        lon_grid_num = (self.MAX_LONGITUDE-self.MIN_LONGITUDE)/self.LON_INTERVAL

        # 判断格网是否是整数
        if  not (lat_grid_num.is_integer() and lon_grid_num.is_integer()):
            raise Exception("please guarantee that the lon/lat_interval can divide 360 or 180")

        # 返回格网个数沿经度方向lon_grid_num, 沿着纬度方向lat_grid_num
        return lon_grid_num,lat_grid_num


    def getGrid(self,shapefile):
        """
        根据h湖北省行政区划图获取栅格，该grid 会标记有经纬度范围，以及ID，根据从下到上，从左到右顺序
        :param shapefile: 湖北省行政区划图矢量，记住，使用之前先做特征融合，只要外部边界
        :return: 返回grid 字典集合，字典格式{ID：（ulx,uly,lrx,lry）}
        """
        # 打开湖北省矢量
        driver = ogr.GetDriverByName('ESRI Shapefile')
        dataSource = driver.Open(shapefile, 0)

        # 记住这里值获取一个layer和第一个feature，
        # 所以正如上面所说，需要将湖北省矢量都聚合了，只要外部边界，只有一个特征，
        # 在envi或者arcgis上面做
        layer = dataSource.GetLayer()
        # 获取投影信息
        self.EPSG = int(layer.GetSpatialRef().GetAuthorityCode("GEOGCS"))
        feature = layer.GetNextFeature()

        # 湖北省矢量边界的轮廓线
        poly_geom = feature.geometry()

        # 湖北省矢量边界的外接矩形
        extent = layer.GetExtent()
        print("湖北省外界矩形",extent)



        # 根据extent 设定 研究区域范围，减小计算量
        self.setUserExtent(extent)

        # print("extent",extent)
        # print(poly_geom)

        gridarr = []
        gridId = 0
        for subgrid in self.user_grid():

            # 用湖北省矢量边界轮廓线生成一个geometry格式的poly
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(subgrid[0],subgrid[1])
            ring.AddPoint(subgrid[2],subgrid[1])
            ring.AddPoint(subgrid[2],subgrid[3])
            ring.AddPoint(subgrid[0],subgrid[3])
            ring.AddPoint(subgrid[0],subgrid[1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            gridId += 1
            if poly.Intersect(poly_geom):
                #记录下这个在湖北省境内的网格在原来矩形网格中的位置
                temp ={"ID":str(gridId),"extent":subgrid}
                gridarr.append(temp)
                temp = None

            ring = None
            poly = None

        return gridarr,self.user_grid()

    def saveGrid(self,shp_path,grid_path):
        # 获取湖北省境内的网格
        mGrid,_ = self.getGrid(shp_path)
        outDriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(grid_path):
            os.remove(grid_path)
        outDataSource = outDriver.CreateDataSource(grid_path)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.EPSG)
        outLayer = outDataSource.CreateLayer(grid_path, srs,geom_type=ogr.wkbPolygon)
        featureDefn = outLayer.GetLayerDefn()
        idField = ogr.FieldDefn("Position", ogr.OFTInteger)
        outLayer.CreateField(idField)
        for subgrid in mGrid:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][1])
            ring.AddPoint(subgrid["extent"][2], subgrid["extent"][1])
            ring.AddPoint(subgrid["extent"][2], subgrid["extent"][3])
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][3])
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            # 生成字段记录下网格在初始湖北省外结矩的位置

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outFeature.SetField("Position", subgrid["ID"])
            outLayer.CreateFeature(outFeature)
            outFeature = None
        outDataSource = None

    def validGrid(self,shp_path,tile_path,output):

        """
        :param shp_path: 湖北省矢量
        :param tile_path: 湖北省Sentinel-2 遥感影像id 矢量
        :param output: 有效npz
        :return: 格王对应哪些遥感影像id
        """
        mGrid,_ = self.getGrid(shp_path)
        driver = ogr.GetDriverByName('ESRI Shapefile')
        dataSource = driver.Open(tile_path, 0)  # 0 means read-only. 1 means writeable.
        if dataSource is None:
            print("Can not open the {}".format(tile_path))
        layer = dataSource.GetLayer()

        if os.path.exists(output):
            os.remove(output)

        gridTilePairs = []
        for subgrid in mGrid:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][1])
            ring.AddPoint(subgrid["extent"][2], subgrid["extent"][1])
            ring.AddPoint(subgrid["extent"][2], subgrid["extent"][3])
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][3])
            ring.AddPoint(subgrid["extent"][0], subgrid["extent"][1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            gridTile = {}
            gridTile["gridId"] = subgrid["ID"]
            gridTile["etRing"] = [subgrid["extent"][0], subgrid["extent"][1],\
                                  subgrid["extent"][2], subgrid["extent"][1],\
                                  subgrid["extent"][2], subgrid["extent"][3],\
                                  subgrid["extent"][0], subgrid["extent"][3],\
                                  subgrid["extent"][0], subgrid["extent"][1]]

            iFeature = 0
            # print(subgrid["ID"])
            tilesintersectwithgrid = []
            for feature in layer:

                geom = feature.geometry()

                if poly.Intersect(geom):
                    tilesintersectwithgrid.append(feature.GetFieldAsString("Name"))
            layer.ResetReading()

            gridTile["valTiles"] = tilesintersectwithgrid # 每个格网对应了几个哨兵-2 的ID
            # print(subgrid["ID"], gridTile)
            gridTilePairs.append(gridTile)


        ID =[gridtile["gridId"] for gridtile in gridTilePairs]
        ValTile = [gridtile["valTiles"] for gridtile in gridTilePairs]
        EtRing =[gridtile["etRing"] for gridtile in gridTilePairs]
        kwargs = {"gridId": ID, "Tile": ValTile,"EtRing":EtRing} # 分别是格网ID 相交的sentinel-2 的ID 以及 格网范围
        # print(ValTile)
        np.savez(output,**kwargs)
























