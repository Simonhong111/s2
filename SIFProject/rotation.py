from osgeo import  gdal,osr,ogr

LL = [['112.8403437215', '29.8271753407887'],\
 ['112.840293165789', '29.8270054058742'],\
 ['112.840224906118', '29.8270227828422'], \
 ['112.840275464746', '29.8271927022589'],
 ['112.8403437215', '29.8271753407887']]

ring = ogr.Geometry(ogr.wkbLinearRing)
for lonlat in LL:
    lon,lat = float(lonlat[0]),float(lonlat[1])
    ring.AddPoint(lon,lat)#将田块的多边形点存入环中

poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)#将环存入多边形中

point = ogr.Geometry(ogr.wkbPoint)
point.AddPoint(112.8403214, 29.8271808)

print(poly.Contains(point))

