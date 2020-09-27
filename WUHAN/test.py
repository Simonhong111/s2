import numpy as np
from matplotlib import pyplot as plt

path  = r'C:\Users\zmhwh\Desktop\Temp\s2latlon.txt'


#read txt method two
f = open(path)
coord = []
for line2 in f:
    coord.append(line2)


coord = [float(c) for c in coord[0].split(" ")[0:-1]]

for i in coord:
    print(i,"end")
x = []
y =[]
print(len(coord))
for i in range(int(len(coord)/2)):

    y.append(coord[2*i])
    x.append(coord[2*i+1])


plt.plot(x,y)
plt.scatter(x,y,s=1,c='r')

points = list(zip(x,y))
print(points)

from shapely.geometry import Polygon
bowtie = Polygon(points)
from shapely.geometry.polygon import LinearRing
lll = LinearRing(points)
print(LinearRing([(1,0), (1,1), (0,0)]))
print("dd",lll.is_ccw)
print("******",bowtie.is_valid)
from shapely.geometry import LineString


clean = bowtie.buffer(0)
print(clean.is_valid)
# help(clean)
# print(clean[0].is_ccw)
exterior = list(clean[0].exterior.coords)
# print(clean)
ex = []
ey = []
for ee in exterior:
    ex.append(float(ee[0]))
    ey.append(float(ee[1]))

ix = []
iy = []
interior = list(clean[1].exterior.coords)
# print("inter",interior)
for ee in interior:
    ix.append(float(ee[0]))
    iy.append(float(ee[1]))

plt.plot(ex,ey,'b*')
plt.plot(ex[0:50],ey[0:50],'r*')
# plt.plot(ix,iy,'r^')
plt.show()