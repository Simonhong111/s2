from osgeo import gdal,osr,ogr
import numpy as np
from matplotlib import pyplot as plt


path = r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif"

raster = gdal.Open(path)
XSize,YSize = raster.RasterXSize,raster.RasterYSize
geo = raster.GetGeoTransform()
proj = osr.SpatialReference()
proj.ImportFromWkt(str(raster.GetProjection()))
print(proj)
src = osr.SpatialReference()
src.ImportFromEPSG(32637)
ct = osr.CoordinateTransformation(proj,src)
p1x,p1y = geo[0],geo[3]
p2x,p2y = geo[0]+XSize*geo[1],geo[3]
p3x,p3y = p2x,geo[3]+YSize*geo[5]
p4x,p4y = p1x,p3y
print(p1x,p1y)
p1x,p1y = ct.TransformPoint(p1y,p1x)[0:2]
p2x,p2y = ct.TransformPoint(p2y,p2x)[0:2]
p3x,p3y = ct.TransformPoint(p3y,p3x)[0:2]
p4x,p4y = ct.TransformPoint(p4y,p4x)[0:2]


ring = ogr.Geometry(ogr.wkbLinearRing)

ring.AddPoint(p1x,p1y)
ring.AddPoint(p2x,p2y)
ring.AddPoint(p3x,p3y)
ring.AddPoint(p4x,p4y)
ring.AddPoint(p1x,p1y)
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

area = poly.Area()/XSize/YSize


print(area)

def Poy2Points(Poly):
    Polystr = str(Poly)
    Polystr = Polystr.split("((")[1]
    Polystr = Polystr.split("))")[0]
    Polystr = Polystr.split(",")
    Points = []
    for p in Polystr:
        if ")" in p:
            p = p[:-1]
        if "(" in p:
            p = p[1:]
        Points.append([float(p.split(" ")[0]),float(p.split(" ")[1])])
    # print("points",Points)
    return np.array(Points)


a = "POLYGON ((119.457430791914 16.6663432028267 0,119.445759332914 16.6914217071463 0,119.455305313832 16.7001576957881 0,119.483468311763 16.7259313599489 0,119.492896420008 16.7431872665541 0,119.432817488393 16.848461115514 0,119.34446854493 16.9853771919713 0,119.325542423192 17.033229034609 0,119.324446371967 17.036000269813 0,119.214870154476 17.0130075252116 0,119.101759348168 16.995770231593 0,119.096871199321 17.0013345156229 0,119.081728231099 16.9877207480642 0,119.020454798759 16.9566702691339 0,118.996887968566 16.9451704709632 0,118.992172077392 16.931365900708 0,119.003944147783 16.8841945506288 0,119.001580460858 16.8508316704706 0,118.994502243349 16.8105669921842 0,118.993319491888 16.7898589365953 0,118.912010927374 16.7530617723507 0,118.883502316397 16.743098075652 0,118.911999685245 16.7053176616213 0,118.962655938833 16.661589059639 0,118.937901486482 16.6167264892547 0,118.932000667469 16.5764616860341 0,118.936693525323 16.48787527416 0,118.933924669009 16.4874770413886 0,118.944940303122 16.4809707006618 0,119.019176818925 16.486707401813 0,119.114621871212 16.4855385471886 0,119.175894995319 16.4832268212432 0,119.245416176464 16.4786136485372 0,119.268980949281 16.4659549752809 0,119.28901223504 16.4625006165515 0,119.305510233561 16.4705514399467 0,119.32082952872 16.4763015734327 0,119.346752827102 16.4751475423318 0,119.393892913692 16.5246111737258 0,119.422182702947 16.603989462924 0,119.44103959695 16.6367753759426 0,119.451647283279 16.6620842838698 0,119.457430791914 16.6663432028267 0,119.457430791914 16.6663432028267 0))"
b = "POLYGON ((119 17 0,120 17 0,120 18 0,119 18 0,119 17 0))"
a = Poy2Points(a)
b = Poy2Points(b)
plt.plot(a[:,0],a[:,1],"r")
plt.plot(b[:,0],b[:,1])
# plt.show()
import csv
path =r"D:\Cornell\EthiopianDrought\CropCSV\Crop\CROP_Index_From2010-2016_V2long.csv"
V = []
W = []

with open(path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row["ID"] != "35763":
            continue
        print("jj")
        for term in ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]:
            V.append(row[term+"2011"])
            W.append(row[term[:-3]+"AREA2011"])

V = np.array(V).astype(np.float)
W = np.array(W).astype(np.float)
print(V)
print(W)
print(np.sum(W*V)/np.sum(W))

print((209.53978*6.39573159+ 198.82709*6.61779082)/ (6.39573159 +6.61779082))