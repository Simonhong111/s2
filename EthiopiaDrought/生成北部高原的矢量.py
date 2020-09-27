from osgeo import ogr,osr
import os



# Create a Polygon from the extent tuple
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(36.4,12.75)
ring.AddPoint(39.6, 12.75)
ring.AddPoint(39.6, 7)
ring.AddPoint(36.4, 7)
ring.AddPoint(36.4,12.75)
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

# Save extent to a new Shapefile
outShapefile = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropSHP\CropLandR.shp"
outDriver = ogr.GetDriverByName("ESRI Shapefile")

# Remove output shapefile if it already exists
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)



# create the spatial reference, WGS84
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)


# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
outLayer = outDataSource.CreateLayer("northhighland", srs,geom_type=ogr.wkbPolygon)

# Add an ID field
idField = ogr.FieldDefn("id", ogr.OFTInteger)
outLayer.CreateField(idField)

# Create the feature and set values
featureDefn = outLayer.GetLayerDefn()
feature = ogr.Feature(featureDefn)
feature.SetGeometry(poly)
feature.SetField("id", 1)
outLayer.CreateFeature(feature)
feature = None

# Save and close DataSource
inDataSource = None
outDataSource = None