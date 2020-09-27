import numpy as np
from matplotlib import pyplot as plt

# path  = r'C:\Users\zmhwh\Desktop\Temp\s2latlon2.txt'

path  = r'C:\Users\zmhwh\Desktop\Temp\validation\s2latlondd32.txt'
#read txt method two
f = open(path)
coord = []
for line2 in f:
    coord.append(line2)
line = coord[0]
line = line.strip().split(" ")

# print(line)
coord = []
for term in line:
    term = term.split(",")
    coord.append(float(term[0]))
    coord.append(float(term[1]))

coord.append(coord[0])
coord.append(coord[1])

# coord = [float(c) for c in coord[0].split(" ")[0:-1]]
#
# for i in coord:
#     print(i,"end")
x = []
y =[]
print(len(coord))
for i in range(int(len(coord)/2)):

    y.append(coord[2*i])
    x.append(coord[2*i+1])
print(x,y)
# print("*",x[-3:],y[-3:])
plt.plot([ 113.8733,113.87273],[30.52,30.520573],'r*')
plt.plot(x,y)
print((113.8733-113.87273)*110000,(30.52-30.520573)*110000)
# i = 4
# plt.plot(x[0:i],y[0:i])
# plt.scatter(x[0:i],y[0:i])
# print(list(zip(x[0:i],y[0:i])))
# plt.plot(x,y)

# plt.scatter(x[195],y[195],s=10,c='r')
# plt.show()
# points = list(zip(x,y))
# # print(points)
#
# from shapely.geometry import Polygon
# bowtie = Polygon(points)
# from shapely.geometry.polygon import LinearRing
# lll = LinearRing(points)
# # print(LinearRing([(1,0), (1,1), (0,0)]))
# print("dd",lll.is_ccw)
# print("******",bowtie.is_valid)
# # from shapely.geometry import LineString
# #
# #
# clean = bowtie.buffer(0)
# print(clean.is_valid)
# # print(clean.exterior.coords)
# # # help(clean)
# # # print(clean[0].is_ccw)
# exterior = list(clean.exterior.coords)
# print(clean)
# ex = []
# ey = []
# for ee in exterior:
#     ex.append(float(ee[0]))
#     ey.append(float(ee[1]))
#
# # ix = []
# # iy = []
# # interior = list(clean[1].exterior.coords)
# # # print("inter",interior)
# # for ee in interior:
# #     ix.append(float(ee[0]))
# #     iy.append(float(ee[1]))
# #
# # plt.plot(ex,ey,'b')
# # plt.plot(ex[0:50],ey[0:50],'r*')
# # plt.plot(ix,iy,'r^')
plt.show()