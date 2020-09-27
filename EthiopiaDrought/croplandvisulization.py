import os
import numpy as np
from osgeo import gdal,osr,ogr
import glob
from matplotlib import pyplot as plt
dirname = r'D:\Cornell\MCD12C1V006Clip'
filenames = glob.glob(os.path.join(dirname,'MCD12C1.A*001.006.*.tif'))
dirname = r'D:\Cornell\MCD12C1v006AggClip'
filenames = glob.glob(os.path.join(dirname,'*MCD12C1.A*001.006.*.tif'))

i = 0
for file in filenames[0:1]:
    i += 1
    fig = plt.figure(i)
    year = os.path.basename(file).split('.')[1][1:5]
    crop = gdal.Open(file).ReadAsArray()
    crop[crop !=12] = 0
    plt.title(year)
    cmap = plt.get_cmap('Greens')

    pim = plt.imshow(crop,cmap=cmap)


plt.show()
