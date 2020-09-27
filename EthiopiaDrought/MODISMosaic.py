from osgeo import gdal,osr,ogr
from datetime import *
from dateutil import rrule
import os
from epdatapro import write_Img
# import datetime
from epdatapro import clipbyshp
import glob
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
def dayofyear(y,m,d):
    """
    :param y:
    :param m:
    :param d:
    :return: the nth day in the year
    """
    days_in_the_year = (date(y, m, d) - date(y, 1, 1)).days + 1
    return str(y)+str(days_in_the_year).zfill(3)

def mosaic(directory):
    start = datetime.strptime("-".join(["2003", '01', "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2020-02-01", "%Y-%m-%d").date()


    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        mdayofyear = dayofyear(dt.year,dt.month,dt.day)
        h21v07 = glob.glob(os.path.join(directory,'*A{}001.h21v07*.hdf'.format(dt.year)))
        h21v08 = glob.glob(os.path.join(directory,'*A{}001.h21v08*.hdf'.format(dt.year)))
        h22v07 = glob.glob(os.path.join(directory, '*A{}001.h22v07*.hdf'.format(dt.year)))
        h22v08 = glob.glob(os.path.join(directory,'*A{}001.h22v08*.hdf'.format(dt.year)))
        assert len(h21v07)*len(h21v08)*len(h22v08)*len(h22v07) >0,"no such file {}".format(dt)
        ph21v07 = gdal.Open(h21v07[0]).GetSubDatasets()[0][0]
        ph21v08 = gdal.Open(h21v08[0]).GetSubDatasets()[0][0]
        ph22v07 = gdal.Open(h22v07[0]).GetSubDatasets()[0][0]
        ph22v08 = gdal.Open(h22v08[0]).GetSubDatasets()[0][0]
        print(ph21v07)
        h21v07evi = gdal.Open(ph21v07).ReadAsArray()
        h21v08evi = gdal.Open(ph21v08).ReadAsArray()
        h22v07evi = gdal.Open(ph22v07).ReadAsArray()
        h22v08evi = gdal.Open(ph22v08).ReadAsArray()

        MosaicImg = np.full(shape=(2400,2400),fill_value=255)
        MosaicImg[0:1200,0:1200] = h21v07evi
        MosaicImg[1200:2400,0:1200] = h21v08evi
        MosaicImg[0:1200, 1200:2400] = h22v07evi
        MosaicImg[1200:2400,1200:2400] = h22v08evi

        preference = gdal.Open(ph21v07)
        proj = osr.SpatialReference()
        proj.ImportFromWkt(str(preference.GetProjection()))
        geotrans = preference.GetGeoTransform()
        path = os.path.join(directory,str(dt.year)+'.tif')

        write_Img(MosaicImg, path, proj, geotrans, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)
        del MosaicImg
        del h22v08evi
        del h21v08evi
        del h21v07evi
        del h22v07evi


def modisclip(modisdirectory,epregionshppath,clippeddirectory):

    modis_files = glob.glob(os.path.join(modisdirectory,"*.tif"))
    # print(chirps_decompress_files)
    for subfile in modis_files:

        input_raster = subfile
        output_raster = os.path.join(clippeddirectory,os.path.basename(subfile))

        # modis = gdal.Open(input_raster)
        # subdateset = modis.GetSubDatasets()[0][0]
        # print(subdateset)
        # gdal.Warp(r"D:\Cornell\EthiopianDrought\Test\2000.02.01all.tif", subdateset, dstSRS='EPSG:4326',
        #           dstNodata=-9999)

        clipbyshp(subfile, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=255)
        print("{} has been processed".format(input_raster))
#
# modisclip(r"D:\Cornell\MCD12C1v006",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\MCD12C1v006AggClip")


# mosaic(r'D:\Cornell\MCD12Q2v006')
# print("I am done")

# datasets = gdal.Open(r"D:\Cornell\MCD12C1v006\MCD12C1.A2003001.006.2018053185458.hdf")
# print(datasets.GetSubDatasets())
# data = datasets.GetSubDatasets()[0][0]
# raster = gdal.Open(data)
# raster_arrlike = raster.ReadAsArray()
# print(raster_arrlike)

def modisclip2(modisdirectory,epregionshppath,clippeddirectory):

    modis_files = glob.glob(os.path.join(modisdirectory,"*.hdf"))
    # print(chirps_decompress_files)
    for subfile in modis_files:

        input_raster = subfile
        output_raster = os.path.join(clippeddirectory,os.path.basename(subfile)+'.tif')

        modis = gdal.Open(input_raster)
        subdateset = modis.GetSubDatasets()[0][0]

        # print(subdateset)


        clipbyshp(subdateset, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=255)
        print("{} has been processed".format(input_raster))

# modisclip2(r"D:\Cornell\MCD12C1v006",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\MCD12C1V006Clip")


# mosaic(r'D:\Cornell\MCD12Q2v006')
# print("I am done")






















