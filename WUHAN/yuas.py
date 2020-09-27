import numpy as np
from matplotlib import pyplot as plt
from osgeo import  gdal,osr,ogr
import pandas as pd
a = [[199980, 3.50004e+06], [409800, 3.50004e+06], [409800, 3.1902e+06], [199980, 3.1902e+06],[199980, 3.50004e+06]]
crs = osr.SpatialReference()
crs.ImportFromEPSG(32650)
wgs = osr.SpatialReference()
wgs.ImportFromEPSG(4326)
ct = osr.CoordinateTransformation(crs,wgs)
points = ct.TransformPoints(a)
print(points)
x1 =[p[1] for p in points]
y1 =[p[0] for p in points]

path  = r'C:\Users\zmhwh\Desktop\Temp\mm4.csv'
data = pd.read_csv(path).to_numpy()
# print(data)
X = data[:,1]
Y = data[:,2]
# print(X,Y)


plt.plot(X,Y)
plt.plot(x1,y1)
plt.show()