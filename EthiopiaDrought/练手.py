from osgeo import osr
import time

ts = time.time()
wgs = osr.SpatialReference()
wgs.ImportFromEPSG(4326)

utm = osr.SpatialReference()
utm.ImportFromEPSG(32637)

for i in range(1000):
    ct = osr.CoordinateTransformation(wgs, utm)
    point = ct.TransformPoint(38, 118)
te = time.time()
print((te-ts)*1000)