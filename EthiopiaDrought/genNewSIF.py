from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
from matplotlib import cm
from matplotlib import pyplot as plt
from epdatapro import write_Img,clipbyshp
def generateNewSIF(sifanomalypath,sifouputdirectory):

    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\SIFAnomaly\sif005_eemd_anomaly_200208RF.nc.tif")

    orgLon = 32.8
    orgLat = 15
    RasterXSize = 310
    RasterYSize = 240
    geotrans = [orgLon,0.05,0.0,orgLat,0.0,-0.05]
    proj = raster.GetProjection()

    sifoutputpath = os.path.join(sifouputdirectory,os.path.basename(sifanomalypath)+".tif")
    eLon = orgLon + RasterXSize*0.05
    eLat = orgLat - RasterYSize*0.05

    orgLon_Col = int((orgLon + 180)*20)
    orgLat_Row = int((90 - orgLat)*20)

    eLon_Col = int((eLon + 180)*20)
    eLat_Row = int((90 - eLat)*20)

    fin = Dataset(sifanomalypath,"r")
    sifanomaly = fin.variables['SIF_740_daily_corr'][:]
    mask = np.where(sifanomaly <= -999)
    sifanomaly[mask] = -9999
    sifanomalydata = sifanomaly[orgLat_Row:eLat_Row,orgLon_Col:eLon_Col]
    write_Img(sifanomalydata, sifoutputpath, proj, geotrans, 310, 240, im_bands=1, dtype=gdal.GDT_Float32)
    fin.close()
    del sifanomalydata
    del sifanomaly

sifanomaly__files = glob.glob(os.path.join(r"D:\Cornell\SIF005_3var\SIF_005", "*.nc"))
# for anomfile in sifanomaly__files:
#     generateNewSIF(anomfile, r"D:\Cornell\NewSIF005TIF")

def sifanomalyClip(sifanomalydirectory,epregionshppath,sifouputdirectory):

    sifanomaly__files = glob.glob(os.path.join(sifanomalydirectory, "*.tif"))
    # print(chirps_decompress_files)
    for anomfile in sifanomaly__files:
        output_raster = os.path.join(sifouputdirectory, os.path.basename(anomfile))
        clipbyshp(anomfile, output_raster, epregionshppath,dstNodata=-9999)
        print("{} has been processed".format(anomfile))

# sifanomalyClip(r"D:\Cornell\NewSIF005TIF",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\NewSIF005Clip")
# sifanomaly__files = glob.glob(os.path.join(r"D:\Cornell\NewSIF005Clip", "*.tif"))
# for file in sifanomaly__files:
#     raster = gdal.Open(file).ReadAsArray()
#     print(raster.max(),raster.min())