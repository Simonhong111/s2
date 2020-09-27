from shapely.geometry import Polygon
from shapely.geometry.polygon import LinearRing
import pandas as pd
from matplotlib import pyplot as plt
path  = r'C:\Users\zmhwh\Desktop\Temp\validation\test.txt'
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
    coord.append([float(term[0]),float(term[1])])


coord.append(coord[0])
print(coord)
data2 = [[p[1],p[0]] for p in coord]
LR = LinearRing(data2)
print("iscc",LR.is_ccw)
# plt.plot(data[:,1],data[:,2])
# plt.show()

coords = [(d[0],d[1]) for d in data2]

bowtie = Polygon(coords)
print(bowtie.is_valid)

