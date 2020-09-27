import os
import glob
import numpy as np
from osgeo import gdal,ogr,osr
"""
主要是根据耕地、水域、湿地等范围，删除掉其中不能用于聚合的东西
"""

def writeLandCoverNpz(lc_path,save_path):
    """
    :param lc_path: 分类数据存放位置
    :return: 写入文件为分类数据影响的地理范围以及该影像的名字
    """
    lc_files = glob.glob(os.path.join(lc_path,"*.tif"))

    lc_name_list = []
    lc_extent_list =[]

    for lc in lc_files:

        raster = gdal.Open(lc)
        GT = raster.GetGeoTransform()
        cols = raster.RasterXSize
        rows = raster.RasterYSize
        lowX = GT[0] + cols*GT[1] + rows*GT[2]
        lowY = GT[3] + cols*GT[4] + rows*GT[5]
        topX = GT[0]
        topY = GT[3]

        lc_name_list.append(lc)
        lc_extent_list.append([topX,topY,lowX,topY,lowX,lowY,topX,lowY,topX,topY])

    kwargs = {"lcName": lc_name_list, "lcExtent": lc_extent_list}
    np.savez(save_path, **kwargs)


def ReadGridNpz(valid_grid_path):
    """
    :param valid_grid_path: 存放了用户自定义的0.05 * 0.05 度的网格以及与其相交的那些Sentinel-2 影像块名称
    :return:返回的为自定义网格对应的ID，对应的Sentinel-2 名称（49SCS），自定义网格的地理范围组成的封闭环（p1->p2->p3->p4->p1）
    """

    gridtilepairs = np.load(valid_grid_path)

    return gridtilepairs["gridId"],gridtilepairs["Tile"],gridtilepairs["EtRing"]

def ReadFilterGridNpz(lc_path):
    """
    :param lc_path: 读取分类数据影像的名称，地理范围组成的环
    :return: 返回分类数据名称，范围组成的环
    """

    mFilterGrid = np.load(lc_path)

    return mFilterGrid["lcName"],mFilterGrid["lcExtent"]

def geo2imgxy(geotransform,xgeo,ygeo):

    trans = geotransform
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([xgeo - trans[0], ygeo - trans[3]])
    return np.linalg.solve(a, b)


def filter(landconver_path,valGridId, valGridTile,valGridEtRing,save_path):
    """
    :param landcoverNames:
    :param landcoverExtent:
    :param valGridId:
    :param valGridEx:
    :param valGridEtRing:
    :return:
    """
    landcover = gdal.Open(landconver_path,0)
    geotransform = landcover.GetGeoTransform()

    lcWidth = landcover.RasterXSize
    lcHeigth = landcover.RasterYSize
    # raster = landcover.ReadAsArray().astype(np.int)
    mGridId =[]
    mGridTile =[]
    mGridRing =[]
    mGridWater=[]
    mGridCropLand=[]
    for idx,mRing in enumerate(valGridEtRing):

        ulx = mRing[0]
        uly = mRing[1]

        lrx = mRing[4]
        lry = mRing[5]

        ulCol,ulRow = geo2imgxy(geotransform,ulx,uly)
        lrCol, lrRow = geo2imgxy(geotransform, lrx, lry)

        width = int(np.ceil(lrCol-ulCol))
        heigth = int(np.ceil(lrRow - ulRow))
        # print(valGridId[idx])
        # print("top",ulCol,ulRow)
        # print("btm",lrCol,lrRow)
        # print(int(ulCol)+ width)
        # print(int(ulRow) + heigth)

        if  (int(ulCol)+ width < lcWidth) and (int(ulRow) + heigth < lcHeigth):

            subImg = landcover.ReadAsArray(int(ulCol), int(ulRow), width, heigth)


            water = np.where(subImg == 60)
            cropland = np.where(subImg == 10)
            waterPercent = (water[0].__len__() * 1.0) / subImg.size
            cropPercent = (cropland[0].__len__() * 1.0) / subImg.size

            if waterPercent < 0.1 and cropPercent > 0.7:
                # print(waterPercent)
                # print(cropPercent)
                # print(subImg.size)
                mGridId.append(valGridId[idx])
                mGridTile.append(valGridTile[idx])
                mGridRing.append(valGridEtRing[idx])
                mGridWater.append(waterPercent)
                mGridCropLand.append(cropPercent)

            del subImg
            del water
            del cropland

        kwargs = {"gridId": mGridId, "gridTile": mGridTile,"EtRing":mGridRing,"waterP":mGridWater,"croplandP":mGridCropLand}
        np.savez(save_path,**kwargs)

def Read_GRID_S2TILE_NPZ(valid_grid_path):
    """
    :param valid_grid_path: 用户自定义网格，所相交Sentinel-2 名称，网格范围组成的环
    :param landcover_path: 分类数据影像名称，范围组成的环
    :param savefilter_path: 过滤出有效的网格，主要是筛选出耕地，无水的网格
    :return:
    """
    valGridId, valGridTile,valGridEtRing = ReadGridNpz(valid_grid_path)

    return valGridId, valGridTile,valGridEtRing


def saveFilterGrid2Shp(filternpz_path,save_shp_path):
    """
    :param filternpz_path: 掩模之后的grid 信息数据 {"gridId": mGridId, "gridTile": mGridTile,"EtRing":mGridRing}
    :param save_shp_path: 将这个矢量存储在什么地方
    :return:
    """

    mfiltergrid = np.load(filternpz_path)
    mGridId = mfiltergrid["gridId"]
    mGridetRing = mfiltergrid["EtRing"]

    outDriver = ogr.GetDriverByName('ESRI Shapefile')

    outDataSource = outDriver.CreateDataSource(save_shp_path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outLayer = outDataSource.CreateLayer(save_shp_path, srs, geom_type=ogr.wkbPolygon)
    featureDefn = outLayer.GetLayerDefn()
    idField = ogr.FieldDefn("Position", ogr.OFTInteger)
    outLayer.CreateField(idField)

    for idx,subgrid in enumerate(mGridetRing):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(subgrid[0], subgrid[1])
        ring.AddPoint(subgrid[2], subgrid[3])
        ring.AddPoint(subgrid[4], subgrid[5])
        ring.AddPoint(subgrid[6], subgrid[7])
        ring.AddPoint(subgrid[8], subgrid[9])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(poly)
        print(mGridId[idx],type(mGridId[idx]))
        outFeature.SetField("Position", mGridId[idx])
        outLayer.CreateFeature(outFeature)
        outFeature = None
    outDataSource = None











