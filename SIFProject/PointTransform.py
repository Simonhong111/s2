from osgeo import gdal,osr,osr
import numpy as np
import os
import glob
from geo2mapxy import *
import xlrd


def fileload(filename = r'D:\HUBEIFIELDDATA\20180711airline\49535149514848535049955750955057_20180711143813.77.xlsx'):
    dataset = []
    workbook = xlrd.open_workbook(filename)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)
    for row in range(table.nrows)[1:]:
        print(table.row_values(row))
        dataset.append(table.row_values(row))
    return dataset,HeadTitle
# data,HeadTitle =fileload(r"H:\temperature_preciption\1998.xlsx")
# print(type(data))

def Wgs84toSinusoidal(dataset,lat,lon):
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)
    target = osr.SpatialReference()
    target.SetProjCS("Sinusoidal")
    target.SetGeogCS("GCS_Unknown", "D_Unknown", "S_Unknown", 6371007.181, 0.0, "Greenwich", 0.0, "Degree")
    target.SetSinusoidal(0, 0, 0)
    ct = osr.CoordinateTransformation(source, target)
    print("lat",lat,"lon",lon)
    geox, geoy = ct.TransformPoint(lat,lon)[:2]
    print("geox", geox, "geoy", geoy)
    PX,PY = geo2imagexy(dataset,geox,geoy)
    return PX,PY

def GetPointByLonLat(dataset,raster,LonLat):

    PX,PY =Wgs84toSinusoidal(dataset,LonLat[0],LonLat[1])
    PX = int(np.floor(PX))
    PY = int(np.floor(PY))
    print(PX,PY)
    return raster[PY,PX]



da = gdal.Open(r"D:\temp\2000058_lsr.tiff")
raster = da.ReadAsArray()
lat = 49+59/60.0+52.97/3600
lon = 100+48/60.0+19.74/3600

value = GetPointByLonLat(da,raster,[32.53,108.53])

print(value)

