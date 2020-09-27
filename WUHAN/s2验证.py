import numpy as np
from matplotlib import pyplot as plt

# path  = r'C:\Users\zmhwh\Desktop\Temp\s2latlon2.txt'

path  = r'C:\Users\zmhwh\Desktop\Temp\validation\ogr.txt'
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

x = []
y =[]
print(len(coord))
for i in range(int(len(coord)/2)):
    y.append(coord[2*i])
    x.append(coord[2*i+1])
print(x,y)
# print("*",x[-3:],y[-3:])
# plt.plot([ 113.8733,113.87273],[30.52,30.520573],'r*')
plt.plot(x,y)
print((113.8733-113.87273)*110000,(30.52-30.520573)*110000)


plt.show()