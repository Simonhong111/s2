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
# print(geometry)

ring = ogr.Geometry(ogr.wkbLinearRing)
# ring.AddPoint(113.83818844073627,31.596390800927402)
# ring.AddPoint(116.04885052164528,31.632004204799575)
# ring.AddPoint( 116.07539603755478,28.836297097147767)
# ring.AddPoint( 113.92627494183455,28.804445245113836)
# ring.AddPoint(113.83818844073627,31.596390800927402)
ring.AddPoint(113.838,31.5964 )
ring.AddPoint(116.049,31.632 )
ring.AddPoint(116.075,28.8363  )
ring.AddPoint(113.926,28.8044 )
ring.AddPoint(113.838,31.5964)

poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

intersection = poly.Intersection(geometry)
# help(intersection)
# help(intersection.GetGeometryRef(0))
print("Boundary",intersection.Boundary())

print("GetGeometryRef",intersection.GetGeometryRef(0))
print("GetGeometryRef",intersection.GetGeometryRef(0).GetX(0))

# print("GetPointCount",intersection.GetPoints())

# print("intersection",intersection)
geometry = extractPointFromFeature(geometry)
X = [p[1] for p in geometry]
Y = [p[0] for p in geometry]

intersection = extractPointFromFeature(intersection)
X1 = [p[1] for p in intersection]
Y1 = [p[0] for p in intersection]

name_dict = {
        'X': X1,
        'Y': Y1
    }

df = pd.DataFrame(name_dict)
df.to_csv(r"C:\Users\zmhwh\Desktop\Temp\mmclip.csv")

# plt.plot(X,Y)
plt.plot(X1,Y1)
print(X1[0:3],Y1[0:3])
plt.plot([ 113.8733,113.87273],[30.52,30.520573],'r*')
print("*",X1[-3:],Y1[-3:])
# plt.show()
