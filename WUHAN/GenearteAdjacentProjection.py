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
ImgPath = r'D:\aria2\hongboshi\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753.SAFE\GRANULE\L2A_T50RKU_A018476_20190105T030703\IMG_DATA\R60m\T50RKU_20190105T030111_B04_60m.jp2'
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
pixel_spacing = 60
p1x,p1y = geo_t[0],geo_t[3]
p2x,p2y = geo_t[0]+60*x_size,geo_t[3]
p3x,p3y = geo_t[0]+60*x_size,geo_t[3]-60*y_size
p4x,p4y = geo_t[0],geo_t[3]-60*y_size
#
WD = [geo_t[0]+i*60+30 for i in range(1830)]
HD = [geo_t[3]-i*60 -30 for i in range(1830)]

Pair = []
for h in HD:
    for w in WD:
        Pair.append([w,h])

Pair= ct.TransformPoints(Pair)
PairX = np.array([p[0] for p in Pair]).reshape(1830,1830)
PairY = np.array([p[1] for p in Pair]).reshape(1830,1830)

PairXPath = r'C:\Users\zmhwh\Desktop\Temp\TestRepX60.0.tif'
PairYPath = r'C:\Users\zmhwh\Desktop\Temp\TestRepY60.0.tif'
write_Img(PairX, PairXPath, src_srs, geo_t,x_size, y_size,im_bands=1, dtype=gdal.GDT_Float32)
write_Img(PairY, PairYPath, src_srs, geo_t,x_size, y_size,im_bands=1, dtype=gdal.GDT_Float32)


# p1x,p1y = ct.TransformPoint(p1x,p1y)[0:2]
# p2x,p2y = ct.TransformPoint(p2x,p2y)[0:2]
# p3x,p3y = ct.TransformPoint(p3x,p3y)[0:2]
# p4x,p4y = ct.TransformPoint(p4x,p4y)[0:2]
#
#
# col,row = 0,0
# px,py = geo_t[0] + col*10,geo_t[3]-10*row
# px,py = ct.TransformPoint(px,py)[0:2]
# print(p2x-p1x)
# print(p4y-p1y)
#
# print('size',x_size,y_size)
# print(col/x_size*(p2x-p1x))
# print(row/y_size*(p4y-p1y))
#
# print("d",p2x-p1x,p3x-p4x)
#
# ctx,cty = p1x + col/x_size*(p3x-p1x),p1y + row/y_size*(p4y-p1y)
#
# print("gdal x {}  -- mine x{}".format(px,ctx),px-ctx)
# print("gdal y {}  -- mine y{}".format(py,cty),py-cty)




# mem_drv = gdal.GetDriverByName("GTiff")
#
# dest = mem_drv.Create(r'C:\Users\zmhwh\Desktop\Temp\tes110.tif', int((p3x - p1x) / pixel_spacing), \
#                       int((p1y - p3y) / pixel_spacing), 1, gdal.GDT_Float32)
#
# new_geo = (p1x, pixel_spacing, geo_t[2], \
#            p1y, geo_t[4], -pixel_spacing)
# # Set the geotransform
# dest.SetGeoTransform(new_geo)
# dest.SetProjection(target_srs.ExportToWkt())
# # Perform the projection/resampling
#
# res = gdal.ReprojectImage(g, dest, \
#                           src_srs.ExportToWkt(), target_srs.ExportToWkt(), \
#                           gdal.GRA_NearestNeighbour)
end_time = time.time()
print("cost times ",end_time-start_time)
# print(res)


