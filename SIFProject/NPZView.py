
import numpy as np



s2 = np.load(r"D:\Sen2Projecton\NPZ\Sen2Agg.npz")

print(s2.files)
['gridId', 'gridEtRing', 'Date', 'mean', 's2path']
s2Id = s2['gridId']
s2Date = s2['Date']
s2mean = s2['mean']

sif = np.load(r"D:\Sen2Projecton\NPZ\2018-5SifAgg.npz")
print(sif.files)
sif_avg = sif["sif757"]
sifid = sif['sifgridId']
sifdate = sif['date']
# print(sif["EtRing"])
# print(set(sifdate))
# print(set(s2Date))
sab = 0
# print(sifid)
# print(s2Id)
# for idx,dt in enumerate(sifdate):
#     for index,tm in enumerate(s2Date):
#         print("*",sifid[idx])
#         print(s2Id[index])

        # if dt == tm and sifid[idx] == s2Id[index]:
        #     print("*")
        #
        #     sab += 1
print(sab)

import os
from osgeo import gdal,osr,ogr
def saveFilterGrid2Shp(filternpz_path,save_shp_path):
    """
    :param filternpz_path: 掩模之后的grid 信息数据 {"gridId": mGridId, "gridTile": mGridTile,"EtRing":mGridRing}
    :param save_shp_path: 将这个矢量存储在什么地方
    :return:
    """

    mfiltergrid = np.load(filternpz_path)
    mGridId = mfiltergrid['sifgridId']
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
        print("***",mGridId[idx])
        outFeature.SetField("Position", str(mGridId[idx]))
        outLayer.CreateFeature(outFeature)
        outFeature = None
    outDataSource = None



saveFilterGrid2Shp(r"D:\Sen2Projecton\NPZ\2018-5SifAgg.npz",r"D:\Sen2Projecton\ShapeFile\SifValid\SifValid.shp")
