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

def chirpsclip(chirpsdirectory,epregionshppath,clippeddirectory,dstNodata=-9999):

    chirps_decompress_files = glob.glob(os.path.join(chirpsdirectory,"*.tif"))
    # print(chirps_decompress_files)
    for tiffile in chirps_decompress_files:
        input_raster = glob.glob(os.path.join(tiffile,"*.tif"))[0]
        output_raster = os.path.join(clippeddirectory,os.path.basename(input_raster))
        clipbyshp(input_raster,output_raster,epregionshppath,dstNodata=dstNodata)
        print("{} has been processed".format(input_raster))


# chirpsclip(r"D:\Cornell\EthiopianDrought\MonthlyRainfall",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\Chirps2")

def calChirpsAverageandStd(chirpsclippeddirectory,outputdirectory,month):


        start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
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

            chirps_file = os.path.join(chirpsclippeddirectory,
                                      "chirps-v2.0."+str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
            print("*", chirps_file)
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

        chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.2010.04.tif")
        geotrans = chirps_reference.GetGeoTransform()
        proj = chirps_reference.GetProjection()
        outputpath = os.path.join(outputdirectory, "chirps_month" + str(month).zfill(2) + ".tif")
        write_Img(data, outputpath, proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)



# for i in range(1,13):
#     calChirpsAverageandStd(r"D:\Cornell\EthiopianDrought\Chirps2",
#                      r"D:\Cornell\EthiopianDrought\AData\Chirps2Pars", i)



def chirpsAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,mm):

    chirpspars_file = os.path.join(chirsparsdirectory,"chirps_month"+str(mm).zfill(2)+".tif")

    if os.path.exists(chirpspars_file):
        print("get the monthly parameters from file {}".format(chirpspars_file))
    else:
        print("can not get the monthly parameters from file {}".format(chirpspars_file))

    start = datetime.strptime("-".join(["2003", str(mm).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    chirps_file = os.path.join(chirpsdirectory, "chirps-v2.0." + str(yy) + "." + str(mm).zfill(2) + ".tif")
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

    anomalyMap[mask2] = (chirps[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap



# yy = "2014"
# mm ="1"
# chirpsdirectory = r"D:\Cornell\EthiopianDrought\Chirps2"
# chirsparsdirectory = r"D:\Cornell\EthiopianDrought\AData\Chirps2Pars"
# chirpsAnomMap = chirpsAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,mm)

def modisclip(modisdirectory,epregionshppath,clippeddirectory):

    modis_files = glob.glob(os.path.join(modisdirectory,"*"))
    # print(chirps_decompress_files)
    for subfile in modis_files:

        input_raster = glob.glob(os.path.join(subfile,"*.hdf"))[0]
        output_raster = os.path.join(clippeddirectory,os.path.basename(subfile)+".tif")
        modis = gdal.Open(input_raster)
        subdateset = modis.GetSubDatasets()[1][0]
        print(subdateset)
        # gdal.Warp(r"D:\Cornell\EthiopianDrought\Test\2000.02.01all.tif", subdateset, dstSRS='EPSG:4326',
        #           dstNodata=-9999)

        clipbyshp(subdateset, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=-3000)
        print("{} has been processed".format(input_raster))

# modisclip(r"D:\Cornell\EthiopianDrought\MOD13C2.006",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia")


def calModisAverageandStd(modisdirectory,outputdirectory,month):
    start = datetime.strptime("-".join(["2003",str(month).zfill(2),"01"]), "%Y-%m-%d").date()
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

        modis_file = os.path.join(modisdirectory,str(dt.year) + "."+str(dt.month).zfill(2)+"."+str(dt.day).zfill(2)+".tif")
        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == -3000)
            maskarr[mask] = -3000
            multidarr[:,:,band_id] = modis
            del modis
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
# for i in range(1,13):
#     calModisAverageandStd(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia",
#                      r"D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyPars", i)

def eviAnomalyMap(evidirectory,eviparsdirectory,yy,mm):

    evi_file = os.path.join(evidirectory,str(yy)+"."+str(mm).zfill(2)+".01"+".tif")
    evipars_file = os.path.join(eviparsdirectory,"evi_month"+str(mm).zfill(2)+".tif")

    evi = gdal.Open(evi_file).ReadAsArray()

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -3000) | (evi == -3000)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -3000) & (evi > -3000)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999
    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def eviAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year) + "." + str(dt.month).zfill(2) +".01" + ".tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year) + "." + str(dt.month).zfill(2) +".01" + ".tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = -3000

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -3000)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average,axisTime

def Pviclip(chirpsdirectory,epregionshppath,clippeddirectory,dstNodata=-9999):

    chirps_decompress_files = glob.glob(os.path.join(chirpsdirectory,"*.tif"))
    # print(chirps_decompress_files)
    for tiffile in chirps_decompress_files:
        input_raster = glob.glob(os.path.join(tiffile,"*.tif"))[0]
        output_raster = os.path.join(clippeddirectory,os.path.basename(input_raster))
        clipbyshp(input_raster,output_raster,epregionshppath,dstNodata=dstNodata)
        print("{} has been processed".format(input_raster))

#
# chirpsclip(r"D:\Cornell\EthiopianDrought\MonthlyRainfall",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\Chirps2")

def calPviAverageandStd(chirpsclippeddirectory,outputdirectory,pvitype):


        start = datetime.strptime("-".join(["2003", str(1).zfill(2), "01"]), "%Y-%m-%d").date()
        stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
        bandNum = 0

        for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

            if pvitype =="all":
                chirps_file = os.path.join(chirpsclippeddirectory,"pvi_{}.tif".format(str(dt.year)))

            else:
                chirps_file = os.path.join(chirpsclippeddirectory,"{}_pvi_{}.tif".format(pvitype,str(dt.year)))

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


def PviAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,pvitype):

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

    anomalyMap[mask2] = (chirps[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap



yy = "2014"
chirpsdirectory = r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018"
chirsparsdirectory = r"D:\Cornell\EthiopianDrought\AData\PviParms"
chirpsAnomMap = PviAnomalyMap(chirpsdirectory,chirsparsdirectory,yy,pvitype="long")
