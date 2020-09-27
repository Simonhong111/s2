import os
import glob
import numpy as np
from osgeo import gdal,ogr,osr

def saveFilterGrid2Shp(sif_path,save_shp_path):
    """
    :param filternpz_path: 掩模之后的grid 信息数据 {"gridId": mGridId, "gridTile": mGridTile,"EtRing":mGridRing}
    :param save_shp_path: 将这个矢量存储在什么地方
    :return:
    """

    mfiltergrid = np.load(sif_path)
    # ['date', 'sif757', 'sif771', 'sifgridId', 'vec_SIF_avg', 'vec_SIF_avg_norm', 'EtRing']
    mGridId = mfiltergrid["sifgridId"]
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

        outFeature.SetField("Position", int(mGridId[idx]))
        outLayer.CreateFeature(outFeature)
        outFeature = None
    outDataSource = None

saveFilterGrid2Shp(r"D:\Satellive\NPZ\2018-5SifAgg.npz",\
                   r"D:\Satellive\SHP\SIFSHP\SIFSHP.shp")