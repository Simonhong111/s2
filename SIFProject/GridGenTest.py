from GridGen import GridDefination
from osgeo import gdal,osr,ogr
import os
import numpy as np


hbshp_path = r"D:\Satellive\SHP\Hubei\Hubei.shp"
tile_path = r"D:\Satellive\SHP\HubeiS2Tile\HubeiS2Tile.shp"



mGridDefination  =   GridDefination()
# _,g= mGridDefination.getGrid(hbshp_path)
# for i in g:
#     print(i)
mGridDefination.validGrid(hbshp_path,tile_path,r"D:\Satellive\NPZ\HB_GRID_S2Tile.npz")
# mGridDefination.saveGrid(hbshp_path,r"D:\Satellive\SHP\HubeiExtentGrid\HubeiExtentGrid.shp")

# for sub in Gridarr:
#     print(sub)
# mGridDefination.validGrid(hbshp_path,tile_path,r"D:\Sen2Projecton\Grid_Tile.npz")
# grid_tile = np.load(r"D:\Sen2Projecton\Grid_Tile.npz")
#
# for idx,item in enumerate(grid_tile["gridId"]):
#     print(item,grid_tile["Tile"][idx],grid_tile["EtRing"][idx])
#
# print(grid_tile.files)






# _,grid =mGridDefination.getGrid(hbshp_path)

# outDriver = ogr.GetDriverByName('ESRI Shapefile')
# if os.path.exists(r"D:\Sen2Projecton\ShapeFile\HBGridAll\HBGridAll.shp"):
#     os.remove(r"D:\Sen2Projecton\ShapeFile\HBGridAll\HBGridAll.shp")
# outDataSource = outDriver.CreateDataSource(r"D:\Sen2Projecton\ShapeFile\HBGridAll\HBGridAll.shp")
# srs = osr.SpatialReference()
# srs.ImportFromEPSG(mGridDefination.EPSG)
# outLayer = outDataSource.CreateLayer(r"D:\Sen2Projecton\ShapeFile\HBGridAll\HBGridAll.shp", srs,geom_type=ogr.wkbPolygon)
# featureDefn = outLayer.GetLayerDefn()
# idField = ogr.FieldDefn("Position", ogr.OFTInteger)
# outLayer.CreateField(idField)
# for idx,subgrid in enumerate(grid):
#     ring = ogr.Geometry(ogr.wkbLinearRing)
#     ring.AddPoint(subgrid[0], subgrid[1])
#     ring.AddPoint(subgrid[2], subgrid[1])
#     ring.AddPoint(subgrid[2], subgrid[3])
#     ring.AddPoint(subgrid[0], subgrid[3])
#     ring.AddPoint(subgrid[0], subgrid[1])
#     poly = ogr.Geometry(ogr.wkbPolygon)
#     poly.AddGeometry(ring)
#     # 生成字段记录下网格在初始湖北省外结矩的位置
#
#     outFeature = ogr.Feature(featureDefn)
#     outFeature.SetGeometry(poly)
#     outFeature.SetField("Position", str(idx+1))
#     outLayer.CreateFeature(outFeature)
#     outFeature = None
# outDataSource = None





# mGridDefination.saveGrid(hbshp_path,r"D:\Sen2Projecton\ShapeFile\HBGrid\HBGrid.shp")



# import matplotlib.pyplot as plt
#
# fig = plt.figure()
#
#
# for subgrid in vhbgrid:
#
#     plt.plot([subgrid["extent"][0], subgrid["extent"][2],subgrid["extent"][2], subgrid["extent"][0],subgrid["extent"][0]],    [subgrid["extent"][1], subgrid["extent"][1],subgrid["extent"][3], subgrid["extent"][3],subgrid["extent"][1]])
# for subgrid in hbgrid:
#     plt.plot([subgrid[0], subgrid[2],subgrid[2], subgrid[0],subgrid[0]],    [subgrid[1], subgrid[1],subgrid[3], subgrid[3],subgrid[1]])
#
# plt.show()


