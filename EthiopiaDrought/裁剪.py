from osgeo import gdal,osr,ogr
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
def clipbyshp(input_raster,output_raster,input_shape, dstNodata=-9999):
    """
    :param input_raster: the raster data being processed later
    :param output_raster: the clipped datas' savepaths
    :param input_shape: the shape defining the extent
    :return: none
    """
    ds = gdal.Warp(output_raster,
                   input_raster,
                   format='GTiff',
                   cutlineDSName=input_shape,  # or any other file format
                   # cutlineDSName=None,
                   # cutlineWhere="FIELD = 'whatever'",
                   # optionally you can filter your cutline (shapefile) based on attribute values
                   cropToCutline=True,
                   dstNodata=dstNodata)  # select the no data value you like
    ds = None
def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

input_raster = r"C:\Users\zmhwh\Downloads\asap_mask_crop_v02.tif"
output_raster = r"D:\Cornell\EthiopianDrought\CropType2015\asap_mask_crop.tif"
input_shape = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"
# clipbyshp(input_raster,output_raster,input_shape, dstNodata=0)

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

# for feature in layer:
#     geom = feature.GetGeometryRef()
#     print (geom,"*********************")

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
plt.imshow(raster,cmap=plt.get_cmap("Greens"))
plt.plot(x,y)
plt.colorbar()
plt.show()

# raster = gdal.Open(r"D:\Cornell\EthiopianDrought\CropType2015\asap_mask_crop.tif")
# geo_t = raster.GetGeoTransform()
# proj = raster.GetProjection()
# W,H = raster.RasterXSize,raster.RasterYSize
# path = r"D:\Cornell\EthiopianDrought\CropType2015\asap_mask_crop_bit.tif"
# data = raster.ReadAsArray()
# data[data >= 120] = 200
# data[data < 120] = 0
# write_Img(data, path, proj, geo_t,W, H,im_bands=1, dtype=gdal.GDT_Float32)