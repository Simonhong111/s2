from osgeo import gdal,osr,ogr
import numpy as np
import os
import time


def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

ImgPath = r'D:\aria2\hongboshi\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753.SAFE\GRANULE\L2A_T50RKU_A018476_20190105T030703\IMG_DATA\R10m\T50RKU_20190105T030111_B04_10m.jp2'
start_time = time.time()

src_srs = osr.SpatialReference()
src_srs.ImportFromEPSG(32650)
target_srs = osr.SpatialReference()
target_srs.ImportFromEPSG(32649)

ct = osr.CoordinateTransformation(src_srs,target_srs)
g = gdal.Open(ImgPath)

geo_t = g.GetGeoTransform()
x_size = g.RasterXSize  # Raster xsize
y_size = g.RasterYSize  # Raster ysize
pixel_spacing = 10
(ulx, uly, ulz) = ct.TransformPoint(geo_t[0], geo_t[3])
(lrx, lry, lrz) = ct.TransformPoint(geo_t[0] + geo_t[1] * x_size, \
                                    geo_t[3] + geo_t[5] * y_size)

mem_drv = gdal.GetDriverByName("GTiff")

dest = mem_drv.Create(r'C:\Users\zmhwh\Desktop\Temp\tesgdal10.tif', int((lrx - ulx) / pixel_spacing), \
                      int((uly - lry) / pixel_spacing), 1, gdal.GDT_Float32)

new_geo = (ulx, pixel_spacing, geo_t[2], \
           uly, geo_t[4], -pixel_spacing)
# Set the geotransform
dest.SetGeoTransform(new_geo)
dest.SetProjection(target_srs.ExportToWkt())
# Perform the projection/resampling

res = gdal.ReprojectImage(g, dest, \
                          src_srs.ExportToWkt(), target_srs.ExportToWkt(), \
                          gdal.GRA_NearestNeighbour)
end_time = time.time()
print("cost times ",end_time-start_time)
print(res)


