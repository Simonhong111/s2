from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob

# class BandComposite(object):
#
#     # 读图像文件
#     def read_img(self, filename):
#         hdf_ds = gdal.Open(filename)  # 打开文件
#         dataset = gdal.Open(hdf_ds.GetSubDatasets()[1][0], gdal.GA_ReadOnly)
#         im_width = dataset.RasterXSize  # 栅格矩阵的列数
#         im_height = dataset.RasterYSize  # 栅格矩阵的行数
#
#
#         im_geotrans = dataset.GetGeoTransform()  # 仿射矩阵
#         im_proj = dataset.GetProjection()  # 地图投影信息
#
#         im_data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 将数据写成数组，对应栅格矩阵
#         mask = (im_data > 65529) | (im_data <0)
#         im_data[mask] = 0
#
#
#         del dataset  # 关闭对象，文件dataset
#         return im_proj, im_geotrans, im_data, im_width, im_height
#
#     # 写文件，以写成tif为例
#     def write_img(self, filename, im_proj, im_geotrans, im_data):
#
#         # 判断栅格数据的数据类型
#         if 'int8' in im_data.dtype.name:
#             # print("int8")
#             datatype = gdal.GDT_Byte
#         elif 'int16' in im_data.dtype.name:
#             # print("int16")
#             datatype = gdal.GDT_UInt16
#         else:
#             # print("int32")
#             datatype = gdal.GDT_Float32
#
#         # 判读数组维数
#         if len(im_data.shape) == 3:
#             im_bands, im_height, im_width = im_data.shape
#         else:
#             im_bands, (im_height, im_width) = 1, im_data.shape
#
#         # 创建文件
#         driver = gdal.GetDriverByName("GTiff")  # 数据类型必须有，因为要计算需要多大内存空间
#         dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)
#
#         dataset.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
#         dataset.SetProjection(im_proj)  # 写入投影
#
#         if im_bands == 1:
#             dataset.GetRasterBand(1).WriteArray(im_data)  # 写入数组数据
#         else:
#             for i in range(im_bands):
#                 dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
#
#         del dataset
#
#
#
#
# files = glob.glob(os.path.join(r"D:\Downloads","*tif"))
#
# for file in files:
#     if os.path.exists(file):
#         os.remove(file)
#
#
# files = glob.glob(os.path.join(r"D:\Downloads","*hdf"))
#
# for file in files:
#     print("**",file)
#     mydata = BandComposite()
#
#     im_proj, im_geotrans, im_data, im_width, im_height = mydata.read_img(file)
#     outfile = file[:-4]+".tif"
#     print("***",outfile)
#     mydata.write_img(outfile, im_proj, im_geotrans, im_data)
#     del mydata
#
# print("************************")
#
# yc_shp_path= r"D:\yichang\yichang.shp"
# driver = ogr.GetDriverByName('ESRI Shapefile')
# dataSource = driver.Open(yc_shp_path, 0)
# layer = dataSource.GetLayer()
# extent = layer.GetExtent()
# londif = 0.1
# latdif = 0.05
# extent =[110.883,110.9667,30.8,30.833]
# print("湖北省外界矩形", extent)
#
#
#
# print("******************************")
# from Sen2Clip import *
#
# scope =[extent[0],extent[3],extent[1],extent[2]]
# scope2 = [extent[0],extent[3],(extent[0]+extent[1])/2,(extent[3]+extent[2])/2]
#
# scope2 = [extent[0]-0.05,extent[3]-0.05,extent[1]-0.05,extent[2]-0.05]
# files = glob.glob(os.path.join(r"D:\Downloads","*tif"))
#
# for file in files:
#     print("generate ",file)
#     outfile = file[:-4]+"_clip.tif"
#     GetRoi(file,outfile,scope)
#
#     outfile = file[:-4] + "_clip2.tif"
#     GetRoi(file, outfile, scope2)
# print("******************************")
#
#
# files = glob.glob(os.path.join(r"D:\Downloads","*_clip.tif"))
#
# Mm = []
# Yy =[]
# i = 0
# for file in files:
#     print("clip",file)
#     raster = gdal.Open(file).ReadAsArray().astype(np.float)
#     mean = np.sum(raster*0.0001)/(len(np.where(raster>0)[0]))
#
#     Yy.append(i)
#     Mm.append(mean)
#     i += 1
#
# print("*************")
#
# files = glob.glob(os.path.join(r"D:\Downloads","*_clip2.tif"))
#
# Mm2 = []
# Yy2 =[]
# i = 0
# for file in files:
#     print("clip2",file)
#     raster = gdal.Open(file).ReadAsArray().astype(np.float)
#     mean = np.sum(raster*0.0001)/(len(np.where(raster)[0]))
#
#     Yy2.append(i)
#     Mm2.append(mean)
#     i += 1
# print("*****************")
# print("Mm",len(Mm))
# print("Mm2",len(Mm2))
# import matplotlib.pyplot as plt
#
# fig = plt.figure()
# print(Mm)
# plt.plot(Yy,Mm)
# plt.plot(Yy2,Mm2,"b")
# Mm3  = np.array(Mm) - np.array(Mm2)
# plt.plot(np.array(Yy2),-Mm3,"r")
# plt.show()

from osgeo import ogr

wkt = "POINT (100 100)"
pt = ogr.CreateGeometryFromWkt(wkt)
bufferDistance = 100
poly = pt.Buffer(bufferDistance,90)

print(str(poly).split(","))
print ("%s buffered by %d is %s" % (pt.ExportToWkt(), bufferDistance, poly.ExportToWkt()))

