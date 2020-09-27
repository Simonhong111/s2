from osgeo import gdal,osr,ogr
import numpy as np
import os
import time
from matplotlib import  pyplot as plt
def write_Img(data, path,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)




    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

# raster = gdal.Open(r"C:\Users\zmhwh\Desktop\beijing.tif").ReadAsArray()
#
# path = r"C:\Users\zmhwh\Desktop\beijing3.tif"
# write_Img(raster, path,9692, 11449,im_bands=1, dtype=gdal.GDT_Byte)

import glob
names = glob.glob(os.path.join(r"D:\MOS\T50RLV","*.tif"))
for name in names:
    print(os.path.basename(name)[-12:])