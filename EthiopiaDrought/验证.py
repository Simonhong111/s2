from osgeo import gdal,osr,ogr
import os
import numpy as np
from matplotlib import pyplot as plt

path = r"D:\Cornell\EthiopianDrought\Chirps2"


# raster = gdal.Open(path).ReadAsArray()



from netCDF4 import Dataset
from osgeo import gdal,osr,ogr

def write_Img(data, path, im_width=3600, im_heigth=1800,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    geotrans = [-180,0.1,0,90,0,-0.1]

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))


    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset


# path =r"C:\Users\zmhwh\Desktop\Temp\ISCCP_HXG_global_radiation_2001_01_01_00.nc"  # 文件路径
# fin = Dataset(path,"r")  # 读入文件
# data = fin.variables['global_radiation'][:] # 读入对应变量
# outpath = r"C:\Users\zmhwh\Desktop\Temp\test.tif" # 生产tif影像，路径
# write_Img(data, outpath, im_width=3600, im_heigth=1800,im_bands=1, dtype=gdal.GDT_Float32)

path = r"C:\Users\zmhwh\Desktop\Temp\ethiopia\Eth_Woreda_2013.shp"

# import geopandas
# myshpfile = geopandas.read_file('myshpfile.shp')
# myshpfile.to_file('myJson.geojson', driver='GeoJSON')

import shapefile
from json import dumps

# read the shapefile
reader = shapefile.Reader(path)
fields = reader.fields[1:]
field_names = [field[0] for field in fields]
buffer = []
for sr in reader.shapeRecords():
    atr = dict(zip(field_names, sr.record))
    geom = sr.shape.__geo_interface__
    buffer.append(dict(type="Feature", \
                       geometry=geom, properties=atr))

# write the GeoJSON file

geojson = open(r"C:\Users\zmhwh\Desktop\Temp\ethiopia\pyshp-demo.json", "w")
geojson.write(dumps({"type": "FeatureCollection", "features": buffer}, indent=2) + "\n")
geojson.close()