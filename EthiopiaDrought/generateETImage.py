import os
import sys
from optparse import OptionParser
import numpy as np
import h5py
import time
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
import time
import glob
from epdatapro import write_Img
from osgeo import gdal,osr,ogr
def generateET(etpath,outputdirectory):

    orgLon = 32
    orgLat = 15
    RasterXSize = 70
    RasterYSize = 50
    geotrans = [orgLon, 0.25, 0.0, orgLat, 0.0, -0.25]

    eLon = orgLon + RasterXSize * 0.25
    eLat = orgLat - RasterYSize * 0.25

    orgLon_Col = int((orgLon + 180) * 4)
    orgLat_Row = int((90 - orgLat) * 4)

    eLon_Col = int((eLon + 180) * 4)
    eLat_Row = int((90 - eLat) * 4)

    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    start = datetime.strptime("1970-01-01", "%Y-%m-%d").date()

    fin = Dataset(etpath, "r")
    obsTime = fin.variables["time"][:]
    et = fin.variables["E"]
    end = timedelta(days=3667) + start



    for id,T in enumerate(obsTime):

        oT = timedelta(days=T) + start
        met = et[id,:,:].T
        data = met[orgLat_Row:eLat_Row,orgLon_Col:eLon_Col]

        path = os.path.join(outputdirectory,
                            "ET.v3.3a" + str(oT.year) + str(oT.month).zfill(2) + ".tif")

        write_Img(data, path, proj, geotrans, 70, 50, im_bands=1, dtype=gdal.GDT_Float32)
        del met
        del data



generateET(r"D:\Cornell\EthiopianDrought\ET\v3.3a\E_1980_2018_GLEAM_v3.3a_MO.nc",r"D:\Cornell\EthiopianDrought\AData\ETMonthV3.3a")

