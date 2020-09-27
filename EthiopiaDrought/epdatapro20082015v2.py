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


def calChirpsAverageandStd(chirpsclippeddirectory,outputdirectory,month):


        start = datetime.strptime("2003-01-31", "%Y-%m-%d").date()
        stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
        bandNum = 0

        for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

            chirps_file = os.path.join(chirpsclippeddirectory,
                                      "chirps-v2.0."+str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")

            if os.path.exists(chirps_file):
                bandNum += 1

        maskarr = np.zeros((228, 299))
        multidarr = np.zeros((228, 299, bandNum))
        band_id = 0

        for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

            if month == "short":
                chirps_file_2 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(2).zfill(2) + ".tif")
                chirps_file_3 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(3).zfill(2) + ".tif")
                chirps_file_4 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(4).zfill(2) + ".tif")
                chirps_file_5 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(5).zfill(2) + ".tif")
            if month == "long":
                chirps_file_2 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(6).zfill(2) + ".tif")
                chirps_file_3 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(7).zfill(2) + ".tif")
                chirps_file_4 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(8).zfill(2) + ".tif")
                chirps_file_5 = os.path.join(chirpsclippeddirectory,
                                             "chirps-v2.0." + str(dt.year) + "." + str(9).zfill(2) + ".tif")

            if os.path.exists(chirps_file_2) and os.path.exists(chirps_file_3) and os.path.exists(chirps_file_4)\
                    and os.path.exists(chirps_file_5):
                chirps2 = gdal.Open(chirps_file_2).ReadAsArray()
                chirps3 = gdal.Open(chirps_file_3).ReadAsArray()
                chirps4 = gdal.Open(chirps_file_4).ReadAsArray()
                chirps5 = gdal.Open(chirps_file_5).ReadAsArray()
                mask = (chirps2 == -9999) | (chirps3 == -9999) | (chirps4 == -9999) | (chirps5 == -9999)
                maskarr[mask] = -9999
                temparr = np.zeros((228, 299, 4))
                temparr[:, :, 0] = chirps2
                temparr[:, :, 1] = chirps3
                temparr[:, :, 2] = chirps4
                temparr[:, :, 3] = chirps5
                print("rainfall temparr",temparr.mean(axis=2).shape)
                multidarr[:, :, band_id] = temparr.mean(axis=2)

                del temparr
                del chirps2
                del chirps3
                del chirps4
                del chirps5
                del mask
                band_id += 1

        averageMatrix = multidarr.mean(axis=2)
        stdMatrix = multidarr.std(axis=2)

        data = np.zeros((228, 299, 3))
        data[:, :, 0] = averageMatrix
        data[:, :, 1] = stdMatrix
        data[:, :, 2] = maskarr

        chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.2010.04.tif")
        geotrans = chirps_reference.GetGeoTransform()
        proj = chirps_reference.GetProjection()
        outputpath = os.path.join(outputdirectory, "chirps_month" + str(month).zfill(2) + ".tif")
        write_Img(data, outputpath, proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)



#
# calChirpsAverageandStd(r"D:\Cornell\EthiopianDrought\Chirps2",
#                      r"D:\Cornell\EthiopianDrought\AData\Chirps2ParsSeanson", "short")



def chirpsAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,mm):


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
    anomalyMap[mask2] = (tempA[mask2] - averageMatrix[mask2])/stdMatrix[mask2]
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



def calModisAverageandStd(modisdirectory,outputdirectory,month):
    start = datetime.strptime("2003-01-01", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,
                                  str(dt.year) + "." + str(dt.month).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr =np.zeros((228,299,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        if month == "short":

            modis_file_2 = os.path.join(modisdirectory,
                                         str(dt.year) + "."+str(2).zfill(2)+"."+str(dt.day).zfill(2)+".tif")
            modis_file_3 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(3).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")
            modis_file_4 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(4).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")
            modis_file_5 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(5).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")

        if month == "long":
            modis_file_2 = os.path.join(modisdirectory,
                                         str(dt.year) + "."+str(6).zfill(2)+"."+str(dt.day).zfill(2)+".tif")
            modis_file_3 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(7).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")
            modis_file_4 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(8).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")
            modis_file_5 = os.path.join(modisdirectory,
                                        str(dt.year) + "." + str(9).zfill(2) + "." + str(dt.day).zfill(2) + ".tif")

        if os.path.exists(modis_file_2) and os.path.exists(modis_file_3) and os.path.exists(modis_file_4) and os.path.exists(modis_file_5):

            modis2 = gdal.Open(modis_file_2).ReadAsArray()*1.0
            modis3 = gdal.Open(modis_file_3).ReadAsArray()*1.0
            modis4 = gdal.Open(modis_file_4).ReadAsArray()*1.0
            modis5 = gdal.Open(modis_file_5).ReadAsArray()*1.0
            mask = (modis2 == -3000) | (modis3 == -3000) | (modis4 == -3000) | (modis5 == -3000)
            maskarr[mask] = -3000
            temparr = np.zeros((228, 299, 4))
            temparr[:, :, 0] = modis2
            temparr[:, :, 1] = modis3
            temparr[:, :, 2] = modis4
            temparr[:, :, 3] = modis5
            print(temparr.shape)
            multidarr[:,:,band_id] = temparr.mean(axis=2)
            del modis2
            del modis3
            del modis4
            del modis5
            del temparr
            del mask
            band_id += 1

    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis= 2)

    data = np.zeros((228,299,3))
    data[:,:,0] = averageMatrix
    data[:,:,1] = stdMatrix
    data[:,:,2] = maskarr

    modis_reference =  gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"evi_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)

# calModisAverageandStd(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia",
#                      r"D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyParsSeason", "short")

def eviAnomalyMap(evidirectory,eviparsdirectory,yy,mm):


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
    anomalyMap[mask2] = (tempA[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap


def calGOSIFAverageandStd(gosifdirectory,outputdirectory,month):
    start = datetime.strptime("2003-01-31", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(gosifdirectory,"{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr =np.zeros((228,299,bandNum))
    band_id = 0


    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        if month == "short":
            gosif_file_2 = os.path.join(gosifdirectory,
                                         "{}GOSIF_Anom{}.tif".format(str(dt.year),str(2).zfill(2)))
            gosif_file_3 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(3).zfill(2)))
            gosif_file_4 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(4).zfill(2)))
            gosif_file_5 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(5).zfill(2)))
        if month == "long":
            gosif_file_2 = os.path.join(gosifdirectory,
                                         "{}GOSIF_Anom{}.tif".format(str(dt.year),str(6).zfill(2)))
            gosif_file_3 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(7).zfill(2)))
            gosif_file_4 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(8).zfill(2)))
            gosif_file_5 = os.path.join(gosifdirectory,
                                        "{}GOSIF_Anom{}.tif".format(str(dt.year), str(9).zfill(2)))
        if os.path.exists(gosif_file_2) and os.path.exists(gosif_file_3) and os.path.exists(gosif_file_4) and os.path.exists(gosif_file_5):

            gosif2 = gdal.Open(gosif_file_2).ReadAsArray()
            gosif3 = gdal.Open(gosif_file_3).ReadAsArray()
            gosif4 = gdal.Open(gosif_file_4).ReadAsArray()
            gosif5 = gdal.Open(gosif_file_5).ReadAsArray()
            mask =(gosif2 == -9999)| (gosif3 == -9999) | (gosif4 == -9999) | (gosif5 == -9999)
            maskarr[mask] = -9999
            temparr = np.zeros((228, 299, 4))
            temparr[:, :, 0] = gosif2
            temparr[:, :, 1] = gosif3
            temparr[:, :, 2] = gosif4
            temparr[:, :, 3] = gosif5
            multidarr[:,:,band_id] = temparr.mean(axis=2)
            del gosif2
            del gosif3
            del gosif4
            del gosif5
            del mask
            del temparr
            band_id += 1

    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis= 2)

    data = np.zeros((228,299,3))
    data[:,:,0] = averageMatrix
    data[:,:,1] = stdMatrix
    data[:,:,2] = maskarr

    modis_reference =  gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"GOSIF_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)

calGOSIFAverageandStd(r"D:\Cornell\EthiopianDrought\AData\GOSIF0317AnomSplit",
                     r"D:\Cornell\EthiopianDrought\AData\GOSIFMonthParsSeason", "long")

def GOSIFAnomalyMap(gosifdirectory,gosifparsdirectory,yy,mm):
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
    anomalyMap[mask2] = (tempA[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap


def fit(y):
    return signal.detrend(y)
def detrend(GOSIFPath,DeGOSIFPath,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(GOSIFPath,
                                   "GOSIF_{}.M{}.tif".format(str(dt.year), str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(GOSIFPath,
                                   "GOSIF_{}.M{}.tif".format(str(dt.year), str(dt.month).zfill(2)))
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            # print("Max: {} Min: {}".format(chirps[chirps > -9999].max(),chirps[chirps > -9999].min()))
            mask = np.where(chirps > 32765)
            maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps

            del chirps
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    multiarrflatten = multidarr.reshape(228*299,bandNum)
    results = map(fit,multiarrflatten)
    deResults = np.array(list(results)).reshape(228,299,bandNum)
    # AVE = []
    for i in range(bandNum):
        deResults[:,:,i][maskarr ==-9999] = -9999
    #     AVE.append(np.mean(multidarr[:,:,i][maskarr > -9999]))
    # return AVE
    modis_reference = gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(DeGOSIFPath, "GOSIF_Anom" + str(month).zfill(2) + ".tif")
    write_Img(deResults, outputpath, proj, geotrans, 299, 228, im_bands=bandNum, dtype=gdal.GDT_Float32)
# Average = []
# for i in range(1,13):
#     res = detrend(r"D:\Cornell\GOSIFV002Clip",
#                      r"D:\Cornell\EthiopianDrought\AData\GOSIF0317Anom",i)
#     Average.append(res)
# print(Average)
# fig = plt.figure()
# for m in range(1,13):
#
#     ax = fig.add_subplot(6, 2, m)
#     ax.set_title("Month_"+str(m).zfill(2))
#     ax.plot([i for i in range(2003,2018)],Average[m-1])
# plt.show()

def detrendSIF(SIFPath,DeSIFPath,month):
    start = datetime.strptime("-".join(["2007", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(SIFPath,
                                   "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(SIFPath,
                                   "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            # print("Max: {} Min: {}".format(chirps[chirps > -9999].max(),chirps[chirps > -9999].min()))
            mask = np.where(chirps == -9999)
            maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps

            del chirps
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    multiarrflatten = multidarr.reshape(228*299,bandNum)
    results = map(fit,multiarrflatten)
    deResults = np.array(list(results)).reshape(228,299,bandNum)
    # AVE = []
    for i in range(bandNum):
        deResults[:,:,i][maskarr ==-9999] = -9999
    #     AVE.append(np.mean(multidarr[:,:,i][maskarr > -9999]))
    # return AVE
    modis_reference = gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_200208.nc.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(DeSIFPath, "NewSIF_Anom" + str(month).zfill(2) + ".tif")
    write_Img(deResults, outputpath, proj, geotrans, 299, 228, im_bands=bandNum, dtype=gdal.GDT_Float32)

# for i in range(1,13):
#     res = detrendSIF(r"D:\Cornell\NewSIF005Clip",
#                      r"D:\Cornell\EthiopianDrought\AData\NewSIF0318Anom2007",i)
#




def calNewSIFAverageandStd(nsifdirectory,outputdirectory,month):
    start = datetime.strptime("2003-01-31", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(nsifdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr =np.zeros((228,299,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        if month == "short":
            nsif_file_2 = os.path.join(nsifdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(2).zfill(2)))
            nsif_file_3 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(3).zfill(2)))
            nsif_file_4 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(4).zfill(2)))
            nsif_file_5 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(5).zfill(2)))
        if month == "long":
            nsif_file_2 = os.path.join(nsifdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(6).zfill(2)))
            nsif_file_3 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(7).zfill(2)))
            nsif_file_4 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(8).zfill(2)))
            nsif_file_5 = os.path.join(nsifdirectory, "{}NewSIF_Anom{}.tif".format(str(dt.year), str(9).zfill(2)))


        if os.path.exists(nsif_file_2) and os.path.exists(nsif_file_3) and os.path.exists(nsif_file_4) and os.path.exists(nsif_file_5):
            nsif2 = gdal.Open(nsif_file_2).ReadAsArray()
            nsif3 = gdal.Open(nsif_file_3).ReadAsArray()
            nsif4 = gdal.Open(nsif_file_4).ReadAsArray()
            nsif5 = gdal.Open(nsif_file_5).ReadAsArray()
            mask = (nsif2 == -9999) | (nsif3 == -9999) | (nsif4 == -9999) | (nsif5 == -9999)
            maskarr[mask] = -9999
            temparr = np.zeros((228, 299, 4))
            temparr[:, :, 0] = nsif2
            temparr[:, :, 1] = nsif3
            temparr[:, :, 2] = nsif4
            temparr[:, :, 3] = nsif5
            multidarr[:,:,band_id] = temparr.mean(axis=2)
            del nsif2
            del nsif3
            del nsif4
            del nsif5
            del mask
            band_id += 1

    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis= 2)

    data = np.zeros((228,299,3))
    data[:,:,0] = averageMatrix
    data[:,:,1] = stdMatrix
    data[:,:,2] = maskarr

    modis_reference =  gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_201512.nc.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"NewSIF_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)

# calNewSIFAverageandStd(r"D:\Cornell\EthiopianDrought\AData\NewSIF0318AnomSplit",
#                      r"D:\Cornell\EthiopianDrought\AData\NewSIFMonthParsSeason", "long")



def NewSIFAnomalyMap(nsifdirectory,nsifparsdirectory,yy,mm):


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

def calPviAverageandStd(chirpsclippeddirectory, outputdirectory, pvitype):

    start = datetime.strptime("-".join(["2003", str(1).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        if pvitype == "all":
            chirps_file = os.path.join(chirpsclippeddirectory, "pvi_{}.tif".format(str(dt.year)))

        else:
            chirps_file = os.path.join(chirpsclippeddirectory, "{}_pvi_{}.tif".format(pvitype, str(dt.year)))

        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        if pvitype == "all":
            chirps_file = os.path.join(chirpsclippeddirectory, "pvi_{}.tif".format(str(dt.year)))
            print("all---", chirps_file)
        else:
            chirps_file = os.path.join(chirpsclippeddirectory, "{}_pvi_{}.tif".format(pvitype, str(dt.year)))
            print(pvitype, "  " + chirps_file)

        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            mask = np.where(chirps == -9999)
            maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps

            del chirps
            del mask
            band_id += 1

    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis=2)

    data = np.zeros((228, 299, 3))
    data[:, :, 0] = averageMatrix
    data[:, :, 1] = stdMatrix
    data[:, :, 2] = maskarr

    chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018\long_pvi_2003.tif")
    geotrans = chirps_reference.GetGeoTransform()
    proj = chirps_reference.GetProjection()
    outputpath = os.path.join(outputdirectory, "pvi_{}".format(pvitype) + ".tif")
    write_Img(data, outputpath, proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)


# calPviAverageandStd(r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018",
#                  r"D:\Cornell\EthiopianDrought\AData\PviParms",pvitype="long")


def PviAnomalyMap(chirpsdirectory, chirsparsdirectory, yy, pvitype):
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

    anomalyMap[mask2] = (chirps[mask2] - averageMatrix[mask2]) / stdMatrix[mask2]

    return anomalyMap
