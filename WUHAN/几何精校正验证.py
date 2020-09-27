from osgeo import gdal,osr,ogr
import numpy as np
import os,glob
from matplotlib import pyplot as plt
# width 12288 heigth 15594
path = r'C:\Users\zmhwh\Desktop\Temp\TY\200322075412_gps.txt'
f = open(path)
coord = []
for line2 in f:
    coord.append(line2)


print(len(coord))
print(len(coord[1].split("\t")))
print(len(coord[10923].split("\t")))

lat = []
lon = []
lon2 = []
lat2 =[]
values = []
for id,item in enumerate(coord):
    if id <2:
        continue
    print("**",id)
    # lon = coord[id].split("\t")[6]
    # lat = coord[id].split("\t")[7]
    print(coord[id][0:7])
    if coord[id][0:7] == "20	3	21":
        print(float(coord[id].split("\t")[6]))
        lon.append(float(coord[id].split("\t")[6]))
        lat.append(float(coord[id].split("\t")[7]))
    if coord[id][0:7] == "20	3	22":
        print(float(coord[id].split("\t")[6]))
        lon2.append(float(coord[id].split("\t")[6])-30)
        lat2.append(float(coord[id].split("\t")[7]))


# print("you",len(lon2))
# plt.plot(lon,lat)
# plt.scatter(lon,lat)
# plt.plot(lon2,lat2)
# plt.scatter(lon2,lat2)
# plt.show()

lineOffset = 2458.924520828017194
sampOffset = 2310.454402137641409
latOffset =   31.002930433653059
longOffset =  114.369840879159639
heightOffset = 42.678974790778796
lineScale = 2816.000000000000000
sampScale = 2688.000000000000000
latScale =    0.237713774635257
longScale =    0.266250235150167
heightScale = 2399.120194648621691

print(114.14758625/2-114.59346044/2)
k = range(-5999,6000,2)
print(len(k))
print( np.s_[0:1])