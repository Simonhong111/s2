from FilteredGridIncudeMore import *
import os
import numpy as np
from osgeo import gdal,osr,ogr
from geo2mapxy import *

valid_grid_path = r"D:\Satellive\NPZ\HB_GRID_S2Tile.npz"  # 这个数据保存了湖北省范围网格，对应的S2遥感影像的ID，以及网格的范围
landcover_path = r"D:\Satellive\HBCrop\HBCrop.tif" # 这个数据保存了湖北省以及周边地区的地物分类栅格地图
Mask_path= r"D:\Satellive\NPZ\MaskedGridMore.npz" # 这个数据保存了用10,220,30,40的信息，种植面积保存的网格信息
Mask_shp_path = r"D:\Satellive\SHP\MaskedGrid\MaskedGridMore.shp"

#
valGridId, valGridTile,valGridEtRing=Read_GRID_S2TILE_NPZ(valid_grid_path)
# #
filter(landcover_path,valGridId, valGridTile,valGridEtRing,Mask_path)
# #
saveFilterGrid2Shp(Mask_path,Mask_shp_path)

#
# maskedgrid = np.load(Mask_path)
# print(maskedgrid.files)
#
# gID = maskedgrid["gridId"]
# gTile = maskedgrid["gridTile"]
# gEtR = maskedgrid["EtRing"]
# gWater = maskedgrid["waterP"]
# gCrop = maskedgrid["croplandP"]
#
#
# print(maskedgrid["gridId"][0],maskedgrid["gridTile"][0],maskedgrid["EtRing"][0])
#
# landcover = gdal.Open(landcover_path, 0)
# geotransform = landcover.GetGeoTransform()


# for idx, mRing in enumerate(gEtR):
#
#
#     from Sen2Clip import *
#
#     GetRoi(landcover_path,r"D:\Satellive\Test\landcoverclip.tif",[mRing[0],mRing[1],mRing[4],mRing[5]])
#     subImg = gdal.Open(r"D:\Satellive\Test\landcoverclip.tif").ReadAsArray()
#
#
#     water = np.where(subImg == 60)
#     cropland = np.where(subImg == 10)
#     waterPercent = (water[0].__len__() * 1.0) / subImg.size
#     cropPercent = (cropland[0].__len__() * 1.0) / subImg.size
#     print("waterPercent",waterPercent,gWater[idx])
#     print("cropPercent",cropPercent,gCrop[idx])
#     print("gWater",gWater[idx])
#     print("gCrop", gCrop[idx])
    # if waterPercent < 0.1 and cropPercent > 0.85:
    #     # print(waterPercent)
    #     # print(cropPercent)
    #     # print(subImg.size)
    #     mGridId.append(valGridId[idx])
    #     mGridTile.append(valGridTile[idx])
    #     mGridRing.append(valGridEtRing[idx])
    #
    # del subImg
    # del water
    # del cropland





