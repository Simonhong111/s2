import os
import glob
from osgeo import gdal
path = r"D:\MOS\T49RGQ"
names = glob.glob(os.path.join(path,"*.tif"))
for name in names:
    raster = gdal.Open(name)
    if raster == None:
        print("no file",name)


