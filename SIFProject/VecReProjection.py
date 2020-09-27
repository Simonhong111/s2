from osgeo import ogr, osr
import os

driver = ogr.GetDriverByName('ESRI Shapefile')

# input SpatialReference
outSpatialRef = osr.SpatialReference()
outSpatialRef = osr.SpatialReference()
outSpatialRef.SetProjCS("Sinusoidal")
outSpatialRef.SetGeogCS("GCS_Unknown","D_Unknown","S_Unknown",6371007.181,0.0,"Greenwich",0.0,"Degree")
outSpatialRef.SetSinusoidal(0,0,0)


# output SpatialReference
inSpatialRef = osr.SpatialReference()
inSpatialRef.ImportFromEPSG(4326)

# create the CoordinateTransformation
coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

# get the input layer
inDataSet = driver.Open(r'D:\UpperYangzeRiverBoundary\Export_Output.shp')
inLayer = inDataSet.GetLayer()

# create the output layer
outputShapefile = r'D:\UpperYangzeRiverBoundary\trans\v2.shp'
if os.path.exists(outputShapefile):
    driver.DeleteDataSource(outputShapefile)
outDataSet = driver.CreateDataSource(outputShapefile)
outLayer = outDataSet.CreateLayer("Sinusoidal_GCS", geom_type=ogr.wkbMultiPolygon)

# add fields
inLayerDefn = inLayer.GetLayerDefn()
for i in range(0, inLayerDefn.GetFieldCount()):
    fieldDefn = inLayerDefn.GetFieldDefn(i)
    outLayer.CreateField(fieldDefn)

# get the output layer's feature definition
outLayerDefn = outLayer.GetLayerDefn()

# loop through the input features
inFeature = inLayer.GetNextFeature()
while inFeature:
    # get the input geometry
    geom = inFeature.GetGeometryRef()
    # reproject the geometry
    geom.Transform(coordTrans)
    # create a new feature
    outFeature = ogr.Feature(outLayerDefn)
    # set the geometry and attribute
    outFeature.SetGeometry(geom)
    for i in range(0, outLayerDefn.GetFieldCount()):
        outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
    # add the feature to the shapefile
    outLayer.CreateFeature(outFeature)
    # dereference the features and get the next input feature
    outFeature = None
    inFeature = inLayer.GetNextFeature()

# Save and close the shapefiles
inDataSet = None
outDataSet = None