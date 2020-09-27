#!/usr/bin/env python

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
def generateESACCI(directory,outputdirectory):
    # Define spatial grid:

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
    start = datetime.strptime("-".join(["1991", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2020-01-01", "%Y-%m-%d").date()



    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        files = glob.glob(
            os.path.join(directory, "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + dt.strftime('%Y%m') + "*fv04.5.nc"))
        nfiles = len(files)
        if nfiles >= 2:

            vec_n = np.zeros((50,70))
            vec_month = np.zeros((50,70))
            for file in files:
                print(file)
                fin = Dataset(file, "r")
                sm_in = fin.variables['sm'][:][0]
                sm_in_selected = sm_in[orgLat_Row:eLat_Row,orgLon_Col:eLon_Col]
                maskarr = np.zeros((50, 70))
                mask = np.where(sm_in_selected > -9999.0)
                maskarr[mask] = 1
                mask = np.where(maskarr == 1)
                vec_n[mask] += 1
                vec_month[mask] = vec_month[mask] + sm_in_selected[mask]
                fin.close()
            # print(vec_n[vec_n < 5])
            vec_n[vec_n <4] = 1
            vec_month[vec_n ==1] = -2
            vec_month = vec_month/vec_n
            print(vec_month.shape)

            path = os.path.join(outputdirectory,"ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-"+str(dt.year)+str(dt.month).zfill(2)+".tif")

            write_Img(vec_month, path, proj, geotrans, 70, 50, im_bands=1, dtype=gdal.GDT_Float32)
            del vec_month
            del vec_n

def generateESACCI2(directory,outputdirectory):
    # Define spatial grid:

    orgLon = 32
    orgLat = 15
    RasterXSize = 70
    RasterYSize = 50
    geotrans = [-180, 0.25, 0.0, 90, 0.0, -0.25]

    eLon = orgLon + RasterXSize * 0.25
    eLat = orgLat - RasterYSize * 0.25

    orgLon_Col = int((orgLon + 180) * 4)
    orgLat_Row = int((90 - orgLat) * 4)

    eLon_Col = int((eLon + 180) * 4)
    eLat_Row = int((90 - eLat) * 4)


    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    start = datetime.strptime("-".join(["1991", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2020-01-01", "%Y-%m-%d").date()

    vec_n = np.zeros((50,70))
    vec_month = np.zeros((50,70))
    files = glob.glob(r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\1992\*")
    data = np.zeros((720,1440,4))
    id = 0
    for file in files:
        fin = Dataset(file, "r")
        lat_in = fin.variables['lat'][:]
        lon_in = fin.variables['lon'][:]
        sm_in = fin.variables['sm'][:][0]
        data[:,:,id] = sm_in
        id += 1






    path = os.path.join(outputdirectory,"test"+".tif")

    write_Img(data, path, proj, geotrans, 1440, 720, im_bands=4, dtype=gdal.GDT_Float32)

#
# generateESACCI2(r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\combined\1992\ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-19920111000000-fv04.5.nc",
#                 r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\ESACCITIF")
# generateESACCI(r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\1992",r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\ESACCITIF")

# directories = glob.glob(r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\combined\*")
# for directory in directories:
#     generateESACCI(directory,
#                    r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\ESACCITIF2")














