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
from scipy import signal
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



def chirpsMap(chirpsdirectory,chirsparsdirectory,yy,mm):


    chirpspars_file = os.path.join(chirsparsdirectory,"chirps_month"+str(mm)+".tif")

    if os.path.exists(chirpspars_file):
        print("get the monthly parameters from file {}".format(chirpspars_file))
    else:
        print("can not get the monthly parameters from file {}".format(chirpspars_file))

    start = datetime.strptime("2003-01-31", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    if mm == "short":
        chirps_file_2 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(2).zfill(2) + ".tif")
        chirps_file_3 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(3).zfill(2) + ".tif")
        chirps_file_4 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(4).zfill(2) + ".tif")
        chirps_file_5 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(5).zfill(2) + ".tif")
    if mm == "long":
        chirps_file_2 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(6).zfill(2) + ".tif")
        chirps_file_3 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(7).zfill(2) + ".tif")
        chirps_file_4 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(8).zfill(2) + ".tif")
        chirps_file_5 = os.path.join(chirpsdirectory,
                                     "chirps-v2.0." + str(yy) + "." + str(9).zfill(2) + ".tif")

    if not (os.path.exists(chirps_file_2) and os.path.exists(chirps_file_3) and os.path.exists(chirps_file_4)\
                    and os.path.exists(chirps_file_5)):
        print("can not find the file {}".format(chirps_file_2))
    chirps2 = gdal.Open(chirps_file_2).ReadAsArray()
    chirps3 = gdal.Open(chirps_file_3).ReadAsArray()
    chirps4 = gdal.Open(chirps_file_4).ReadAsArray()
    chirps5 = gdal.Open(chirps_file_5).ReadAsArray()
    chirpspars = gdal.Open(chirpspars_file)

    averageMatrix = chirpspars.GetRasterBand(1).ReadAsArray()
    stdMatrix = chirpspars.GetRasterBand(2).ReadAsArray()
    maskarr = chirpspars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (chirps2 == -9999) | (chirps3 == -9999) | (chirps4 == -9999) | (chirps5 == -9999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (chirps2 > -9999) & (chirps3 > -9999) & (chirps4 > -9999) & (chirps5 > -9999)
    mask2 = np.where(mask2)

    anomalyMap = np.zeros(shape=chirps2.shape)
    anomalyMap[mask] = -9999

    temparr = np.zeros((228, 299, 4))
    temparr[:, :, 0] = chirps2
    temparr[:, :, 1] = chirps3
    temparr[:, :, 2] = chirps4
    temparr[:, :, 3] = chirps5
    tempA = temparr.mean(axis=2)
    anomalyMap[mask2] = tempA[mask2]
    del temparr
    return anomalyMap



# yy = "2014"
# mm ="short"
# chirpsdirectory = r"D:\Cornell\EthiopianDrought\Chirps2"
# chirsparsdirectory = r"D:\Cornell\EthiopianDrought\AData\Chirps2ParsSeanson"
# chirpsAnomMap = chirpsAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,mm)
# cax = plt.imshow(chirpsAnomMap,cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
# cbar = plt.colorbar(cax,fraction=0.036,pad=0.04)
# plt.show()


def eviMap(evidirectory,eviparsdirectory,yy,mm):


    evipars_file = os.path.join(eviparsdirectory,"evi_month"+str(mm).zfill(2)+".tif")

    if mm == "short":
        evi_file_2 = os.path.join(evidirectory,
                                    str(yy) + "." + str(2).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_3 = os.path.join(evidirectory,
                                  str(yy) + "." + str(3).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_4 = os.path.join(evidirectory,
                                  str(yy) + "." + str(4).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_5 = os.path.join(evidirectory,
                                  str(yy) + "." + str(5).zfill(2) + "." + str(1).zfill(2) + ".tif")
    if mm == "long":
        evi_file_2 = os.path.join(evidirectory,
                                    str(yy) + "." + str(6).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_3 = os.path.join(evidirectory,
                                  str(yy) + "." + str(7).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_4 = os.path.join(evidirectory,
                                  str(yy) + "." + str(8).zfill(2) + "." + str(1).zfill(2) + ".tif")
        evi_file_5 = os.path.join(evidirectory,
                                  str(yy) + "." + str(9).zfill(2) + "." + str(1).zfill(2) + ".tif")




    evi2 = gdal.Open(evi_file_2).ReadAsArray() * 1.0
    evi3 = gdal.Open(evi_file_3).ReadAsArray() * 1.0
    evi4 = gdal.Open(evi_file_4).ReadAsArray() * 1.0
    evi5 = gdal.Open(evi_file_5).ReadAsArray() * 1.0

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -3000) | (evi2 == -3000) | (evi3 == -3000) | (evi4 == -3000) |(evi5 == -3000)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -3000) & (evi2 > -3000) & (evi3 > -3000) & (evi4 > -3000) & (evi5 > -3000)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi2.shape)
    anomalyMap[mask] = -9999
    temparr = np.zeros((228, 299, 4))
    temparr[:, :, 0] = evi2
    temparr[:, :, 1] = evi3
    temparr[:, :, 2] = evi4
    temparr[:, :, 3] = evi5
    tempA = temparr.mean(axis=2)
    anomalyMap[mask2] = tempA[mask2]

    return anomalyMap




def GOSIFMap(gosifdirectory,gosifparsdirectory,yy,mm):
    if mm == "short":
        gosif_file_2 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(2).zfill(2)))
        gosif_file_3 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(3).zfill(2)))
        gosif_file_4 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(4).zfill(2)))
        gosif_file_5 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(5).zfill(2)))
    if mm == "long":
        gosif_file_2 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(6).zfill(2)))
        gosif_file_3 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(7).zfill(2)))
        gosif_file_4 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(8).zfill(2)))
        gosif_file_5 = os.path.join(gosifdirectory,
                                    "{}GOSIF_Anom{}.tif".format(str(yy), str(9).zfill(2)))

    gosifpars_file = os.path.join(gosifparsdirectory,"GOSIF_month"+str(mm).zfill(2)+".tif")

    gosif2 = gdal.Open(gosif_file_2).ReadAsArray()
    gosif3 = gdal.Open(gosif_file_3).ReadAsArray()
    gosif4 = gdal.Open(gosif_file_4).ReadAsArray()
    gosif5 = gdal.Open(gosif_file_5).ReadAsArray()


    gosifpars = gdal.Open(gosifpars_file)
    averageMatrix = gosifpars.GetRasterBand(1).ReadAsArray()
    stdMatrix = gosifpars.GetRasterBand(2).ReadAsArray()
    maskarr = gosifpars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (gosif2 == -9999) | (gosif2 == -9999) | (gosif2 == -9999) | (gosif2 == -9999)\
    | (gosif2 >32765) | (gosif3 >32765) | (gosif4 >32765) | (gosif5 >32765)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (gosif2 > -9999) & (gosif3 > -9999) & (gosif4 > -9999) & (gosif5 > -9999)\
    & (gosif2 < 32766) & (gosif3 < 32766) & (gosif4 < 32766) & (gosif5 < 32766)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=gosif2.shape)
    anomalyMap[mask] = -9999
    temparr = np.zeros((228, 299, 4))
    temparr[:, :, 0] = gosif2
    temparr[:, :, 1] = gosif3
    temparr[:, :, 2] = gosif4
    temparr[:, :, 3] = gosif5
    tempA = temparr.mean(axis=2)
    anomalyMap[mask2] = tempA[mask2]

    return anomalyMap



def NewSIFMap(nsifdirectory,nsifparsdirectory,yy,mm):


    if mm == "short":
        nsif_file_2 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(2).zfill(2)))
        nsif_file_3 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(3).zfill(2)))
        nsif_file_4 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(4).zfill(2)))
        nsif_file_5 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(5).zfill(2)))
    if mm == "long":
        nsif_file_2 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(6).zfill(2)))
        nsif_file_3 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(7).zfill(2)))
        nsif_file_4 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(8).zfill(2)))
        nsif_file_5 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(yy), str(9).zfill(2)))
    nsifpars_file = os.path.join(nsifparsdirectory,"NewSIF_month"+str(mm).zfill(2)+".tif")

    nsif2 = gdal.Open(nsif_file_2).ReadAsArray()
    nsif3 = gdal.Open(nsif_file_3).ReadAsArray()
    nsif4 = gdal.Open(nsif_file_4).ReadAsArray()
    nsif5 = gdal.Open(nsif_file_5).ReadAsArray()

    nsifpars = gdal.Open(nsifpars_file)
    averageMatrix = nsifpars.GetRasterBand(1).ReadAsArray()
    stdMatrix = nsifpars.GetRasterBand(2).ReadAsArray()
    maskarr = nsifpars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (nsif2 == -9999) | (nsif3 == -9999) | (nsif4 == -9999) | (nsif5 == -9999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (nsif2 > -9999) & (nsif3 > -9999) & (nsif4 > -9999) & (nsif5 > -9999)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=nsif2.shape)
    anomalyMap[mask] = -9999
    temparr = np.zeros((228, 299, 4))
    temparr[:, :, 0] = nsif2
    temparr[:, :, 1] = nsif3
    temparr[:, :, 2] = nsif4
    temparr[:, :, 3] = nsif5
    tempA = temparr.mean(axis=2)
    anomalyMap[mask2] = (tempA[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap


def PviMap(chirpsdirectory, chirsparsdirectory, yy, pvitype):
    if pvitype == "all":
        chirpspars_file = os.path.join(chirsparsdirectory, "pvi_all.tif")
    else:
        chirpspars_file = os.path.join(chirsparsdirectory, "pvi_{}.tif".format(pvitype))

    if os.path.exists(chirpspars_file):
        print("get the monthly parameters from file {}".format(chirpspars_file))
    else:
        print("can not get the monthly parameters from file {}".format(chirpspars_file))

    start = datetime.strptime("-".join(["2003", str(1).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    if pvitype == "all":
        chirps_file = os.path.join(chirpsdirectory, "pvi_{}.tif".format(str(yy)))
        print("all---", chirps_file)
    else:
        chirps_file = os.path.join(chirpsdirectory, "{}_pvi_{}.tif".format(pvitype, str(yy)))
        print(pvitype, "  " + chirps_file)

    if not os.path.exists(chirps_file):
        print("can not find the file {}".format(chirps_file))

    chirps = gdal.Open(chirps_file).ReadAsArray()
    chirpspars = gdal.Open(chirpspars_file)
    averageMatrix = chirpspars.GetRasterBand(1).ReadAsArray()
    stdMatrix = chirpspars.GetRasterBand(2).ReadAsArray()
    maskarr = chirpspars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (chirps == -9999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (chirps > -9999)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=chirps.shape)
    anomalyMap[mask] = -9999

    anomalyMap[mask2] = chirps[mask2]

    return anomalyMap
