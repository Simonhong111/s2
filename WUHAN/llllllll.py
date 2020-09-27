from osgeo import gdal,osr,ogr
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
def extractPointFromFeature(feature):
    feature = str(feature)
    polyNum = len(feature.split("))"))
    print("feature number is ",polyNum)
    feature = feature.split("((")[1].split("))")[0].split(",")
    points = [[float(f.split(" ")[1]),float(f.split(" ")[0])] for f in feature]
    return points
geojson = ogr.Open(r'C:\Users\zmhwh\Desktop\Temp\wuhan-new.json')
layer = geojson.GetLayer()
feature = layer.GetFeature(0)
geometry = feature.GetGeometryRef()
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(113.838,31.5964 )
ring.AddPoint(116.049,31.632 )
ring.AddPoint(116.075,28.8363 )
ring.AddPoint(113.926,28.8044 )
ring.AddPoint(113.838,31.5964)
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

intersection = poly.Intersection(geometry)

print("Boundary",intersection.Boundary())

print("GetGeometryRef",intersection.GetGeometryRef(0))
print("GetGeometryRef",intersection.GetGeometryRef(0).GetX(0))

# print("GetPointCount",intersection.GetPoints())

boundary = intersection.Boundary()
boundaryNum = boundary.GetPointCount()
bX =[]
bY =[]
for i in range(boundaryNum):
    bX.append(boundary.GetX(i))
    bY.append(boundary.GetY(i))

geom = intersection.GetGeometryRef(0)
geomNum = geom.GetPointCount()
gX =[]
gY =[]
for i in range(geomNum):
    gX.append(boundary.GetX(i))
    gY.append(boundary.GetY(i))

geometry = extractPointFromFeature(geometry)
X = [p[1] for p in geometry]
Y = [p[0] for p in geometry]

intersection = extractPointFromFeature(intersection)
X1 = [p[1] for p in intersection]
Y1 = [p[0] for p in intersection]


plt.plot(bX,bY,"r*")
plt.plot(gX,gY,'g*')
plt.plot(X1,Y1)
plt.show()

gX = np.array(gX)
gY = np.array(gY)
bX = np.array(bX)
bY = np.array(bY)
print("gx-bx",gX-bX)
print("gy - by",gY-bY)

X1 =np.array(X1)
Y1 =np.array(Y1)
print("gx-x1",gX-X1)
print("gy-y1",gY-Y1)
