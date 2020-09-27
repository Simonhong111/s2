from shapely.geometry import Polygon
from shapely.geometry.polygon import LinearRing
import pandas as pd
from matplotlib import pyplot as plt
path  = r'C:\Users\zmhwh\Desktop\Temp\mmclip.csv'
data = pd.read_csv(path).to_numpy()
# print(data)
X = data[:,1]
Y = data[:,2]
# print(X,Y)
coord = []
for i,term in enumerate(X):
    coord.append([term,Y[i]])

PP = [113.867982511734,
30.6952740476922,
113.899179492638,
29.7058729412524,
115.033286585814,
29.7275896486434,
115.013476249921,
30.7178664204507,
113.867982511734,
30.6952740476922]

X2= [PP[2*i] for i in range(int(len(PP)/2))]
Y2= [PP[2*i+1] for i in range(int(len(PP)/2))]

data2 = [[p[1],p[0]] for p in coord]
LR = LinearRing(data2)
print("iscc",LR.is_ccw)
plt.plot(X,Y)
plt.plot(X2,Y2)
plt.show()

coords = [(d[0],d[1]) for d in data2]

bowtie = Polygon(coords)
print(bowtie.is_valid)

