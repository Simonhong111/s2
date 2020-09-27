from osgeo import gdal,osr,ogr
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

input_shape = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"
rainpath = r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip.tif"
rain = gdal.Open(rainpath)
geo_t = rain.GetGeoTransform()
print(geo_t)
daShapefile = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()

geo = str(geo).split("((")[1].split("))")[0].split(",")
x = []
y = []
for term in geo:
    x.append(float(term.split(" ")[0]))
    y.append(float(term.split(" ")[1]))

x = np.array(x)
y = np.array(y)
x = (x - geo_t[0]) / geo_t[1]
y = (y - geo_t[3]) / geo_t[5]

raster = gdal.Open(rainpath).ReadAsArray()*1.0

raster[raster == 255] = 255
raster[raster != 255] = np.nan
plt.imshow(raster)
plt.plot(x,y)
plt.show()

