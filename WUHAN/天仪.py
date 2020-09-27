import os,glob
import numpy as np
from osgeo import gdal
from matplotlib import pyplot as plt


path = r"D:\BaiduNetdiskDownload\天仪研究院遥感样例\原始影像\0807\190804160133_10006_b_UTC.raw"
size = os.path.getsize(path)
if np.mod(size,32000) == 0:
    column = 16000
if np.mod(size,24480) == 0:
    column = 12240
row = int(size/(column*2))
print(column,row)


# #!/usr/bin/env python2
# """A simple (slightly broken) implementation of writing a VRT file from Python, using GDAL.
# Adapted from sample code at: http://gis.stackexchange.com/questions/44003/python-equivalent-of-gdalbuildvrt#answer-44048
# """
# ## Options ##
# x_size, y_size = column, row
# source_path = r'D:\BaiduNetdiskDownload\0717\190714031139_15000_UTC.raw'
# source_band = 1
# x_source_size, y_source_size = column, row
# x_offset, y_offset = 0, 0
# dest_x_offset, dest_y_offset = 0, 0
# x_dest_size, y_dest_size = column, row
#
# drv = gdal.GetDriverByName("VRT")
# vrt = drv.Create(r"D:\BaiduNetdiskDownload\0717\190714031139_15000_UTC.vrt", x_size, y_size, 0)
#
# vrt.AddBand(gdal.GDT_UInt16)
# band = vrt.GetRasterBand(1)
#
#
# simple_source = '<SourceFilename relativeToVRT="1">%s</SourceFilename>' % source_path + \
#     '<SourceBand>%i</SourceBand>' % source_band + \
#     '<SourceProperties RasterXSize="%i" RasterYSize="%i" DataType="Real"/>' % (x_source_size, y_source_size) + \
#     '<SrcRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (x_offset, y_offset, x_source_size, y_source_size) + \
#     '<DstRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (dest_x_offset, dest_y_offset, x_dest_size, y_dest_size)
#
# band.SetMetadataItem("SimpleSource", simple_source)
#
# # Changed from an integer to a string, since only strings are allowed in `SetMetadataItem`.
# band.SetMetadataItem("NoDataValue", '-9999')


# gdal.Translate(r"D:\BaiduNetdiskDownload\0717\190714031139_15000_UTC.tif",r"D:\BaiduNetdiskDownload\0717\190714031139_15000_UTC.vrt")

# raster = gdal.Open(r"D:\BaiduNetdiskDownload\天仪研究院遥感样例\原始影像\0807\190804160133_10006_b_UTC.raw")
# data = raster.ReadAsArray(0,0,5000,5000)
# plt.imshow(data)

# lon = [-77.150,-77.150,-77.033,-76.967			]
# lat =[36.417,36.417,36.367,36.333]
# plt.show()
