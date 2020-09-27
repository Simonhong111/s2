import numpy as np
import os
from osgeo import gdal,ogr,osr
import matplotlib.pyplot as plt

def ReadVertexFromTileFile(shp_path):

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shp_path, 0)
    if dataSource is None:
        print('Could not open {}'.format(shp_path))
        return
    else:
        print('Opened {}'.format(shp_path))
        layer = dataSource.GetLayer()
        featureCount = layer.GetFeatureCount()
        print("Number of features in {}: {}".format(os.path.basename(shp_path), featureCount))
    GridInfo = []
    for feature in layer:
        tileinfo_dict = {}
        tileinfo_dict["TileIDName"] = feature.GetFieldAsString("Name")
        print(feature.GetFieldAsString("Name"))
        geom = feature.GetGeometryRef()
        # print(help(geom))
        print()
        tileinfo_dict["CenterLon"] = float(str(geom.Centroid())[7:-1].split(" ")[0])
        tileinfo_dict["CenterLat"] = float(str(geom.Centroid())[7:-1].split(" ")[1])
        tileinfo_dict["TileLon"],tileinfo_dict["TileLat"] = GetPointFromGeom(geom)
        # print(GetPointFromGeom(geom))

        GridInfo.append(tileinfo_dict)

        # print((str(geom)[10:-2]).split(","))
    layer.ResetReading()
    return GridInfo
def GetPointFromGeom(geometry):
    geom = str(geometry)[10:-2].split(",")
    pointLon = []
    pointLat = []
    for g in geom:

        pointLon.append(float(g.split(" ")[0]))
        pointLat.append(float(g.split(" ")[1]))
    return pointLon,pointLat


def GetShapeFileField(shp_path):

    dataSource = ogr.Open(shp_path)
    daLayer = dataSource.GetLayer(0)
    layerDefinition = daLayer.GetLayerDefn()

    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())

def ReadVertexFromShapeFile(shp_path):

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shp_path, 0)
    if dataSource is None:
        print('Could not open {}'.format(shp_path))
        return
    else:
        print('Opened {}'.format(shp_path))
        layer = dataSource.GetLayer()
        featureCount = layer.GetFeatureCount()
        print("Number of features in {}: {}".format(os.path.basename(shp_path), featureCount))
    GridInfo = []
    for feature in layer:

        tileinfo_dict = {}
        # tileinfo_dict["TileIDName"] = feature.GetFieldAsString("Name")
        # print(feature.GetFieldAsString("Name"))
        geom = feature.GetGeometryRef()
        tileinfo_dict["TileLon"],tileinfo_dict["TileLat"] = GetPointFromGeom(geom)
        # print(geom)
        # print(GetPointFromGeom(geom))
        GridInfo.append(tileinfo_dict)
        #
        # print((str(geom)[10:-2]).split(","))
    layer.ResetReading()
    return GridInfo




# hbs2path = r"D:\Sen2Projecton\HBS2TILE\HBS2.shp"
# hbshppath = r"D:\Sen2Projecton\HBSHP\Hubei.shp"
#
# # GetShapeFileField(hbs2path)
# mGridInfo = ReadVertexFromTileFile(hbs2path)
# mShapeFile = ReadVertexFromShapeFile(hbshppath)
#
# fig = plt.figure()
#
# for mgrid in mGridInfo:
#     lon = mgrid["TileLon"]
#     lat = mgrid["TileLat"]
#     centerLon = mgrid["CenterLon"]
#     centerLat = mgrid["CenterLat"]
#
#     plt.plot(lon,lat)
#     plt.scatter(centerLon,centerLat)
#     plt.text(centerLon,centerLat,mgrid["TileIDName"])
# lon = mShapeFile[0]["TileLon"]
# lat = mShapeFile[0]["TileLat"]
#
# plt.plot(lon,lat)
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")
# plt.show()
