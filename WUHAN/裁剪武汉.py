from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob
import time
# from matplotlib import pyplot as plt
# x=[120.97721121345067,120.927801655897,120.95502733046739,121.0246040543695,120.97721121345067]
# y = [46.0887206357723,46.09729168147039,46.122500639405935,46.114937952025265,46.0887206357723]
# x1=[114.79992,115.439101,115.300697,114.664411]
# y1 = [24.424296,24.297468,23.717025,23.843602]
#
# plt.plot(x,y)
# plt.plot(x1,y1,'r')
# plt.show()

rs = np.array([3,6,10,11,14,19,40])
sa =np.array([22,13,17,43,24,23,85])
zb = rs/sa*100
print(zb)
