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

#
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


def chirpsAveMap(chirpsdirectory,yy, month):

# average of a specifical image
    chirps_file = os.path.join(chirpsdirectory,
                               "chirps-v2.0." + str(yy) + "." + str(month).zfill(2) + ".tif")

    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip',"MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file) and os.path.exists(cropland):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -9999) & (crop == 12)
        ave = np.mean(chirps[mask])
        return ave, month



def calChirpsAllAve(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]

        if os.path.exists(chirps_file) and os.path.exists(cropland):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop !=12] = -9999
            del crop
            del chirps
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    ave = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std


def chirpsAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "chirps-v2.0." + str(dt.year) + "." + str(dt.month).zfill(2) + ".tif")

        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]


        # print("*", chirps_file)
        if os.path.exists(chirps_file) and os.path.exists(cropland):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = -9999
            del crop
            del chirps
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)

    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average, axisTime


def generateSIFAnomaly(sifanomalypath,sifouputdirectory):

    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.2010.04.tif")
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
    sifanomaly = fin.variables['SIF_740_daily_corr_anomaly'][:]
    mask = np.where(sifanomaly <= -9999)
    sifanomaly[mask] = -9999
    sifanomalydata = sifanomaly[orgLat_Row:eLat_Row,orgLon_Col:eLon_Col]
    write_Img(sifanomalydata, sifoutputpath, proj, geotrans, 310, 240, im_bands=1, dtype=gdal.GDT_Float32)
    fin.close()
    del sifanomalydata
    del sifanomaly

# sifanomaly__files = glob.glob(os.path.join(r"D:\Cornell\EthiopianDrought\Detrend_RFanom", "*.nc"))
# for anomfile in sifanomaly__files:
#     generateSIFAnomaly(anomfile, r"D:\Cornell\EthiopianDrought\SIFAnomaly")

def sifanomalyClip(sifanomalydirectory,epregionshppath,sifouputdirectory):

    sifanomaly__files = glob.glob(os.path.join(sifanomalydirectory, "*.tif"))
    # print(chirps_decompress_files)
    for anomfile in sifanomaly__files:
        output_raster = os.path.join(sifouputdirectory, os.path.basename(anomfile))
        clipbyshp(anomfile, output_raster, epregionshppath)
        print("{} has been processed".format(anomfile))

# sifanomalyClip(r"D:\Cornell\EthiopianDrought\SIFAnomaly",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\clipedSIf")



def calSIFAverageandStd(chirpsclippeddirectory, outputdirectory, month):

    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsclippeddirectory,
                                   "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsclippeddirectory,
                                    "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")

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
    # print(stdMatrix)
    chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\clipedSIf\sif005_eemd_anomaly_200208RF.nc.tif")
    geotrans = chirps_reference.GetGeoTransform()
    proj = chirps_reference.GetProjection()
    outputpath = os.path.join(outputdirectory, "sif_month" + str(month).zfill(2) + ".tif")
    write_Img(data, outputpath, proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)


# for i in range(1,13):
#     calSIFAverageandStd(r"D:\Cornell\EthiopianDrought\clipedSIf",
#                      r"D:\Cornell\EthiopianDrought\AData\SIFPars", i)

def sifAnomalyMap2(chirpsdirectory,chirsparsdirectory,yy,mm):

    chirpspars_file = os.path.join(chirsparsdirectory,"sif_month"+str(mm).zfill(2)+".tif")
    if not os.path.exists(chirpspars_file):
        print("can not find the monthly parameters file {}".format(chirpspars_file))

    start = datetime.strptime("-".join(["2003", str(mm).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    chirps_file = os.path.join(chirpsdirectory,"sif005_eemd_anomaly_" + str(yy) + str(mm).zfill(2) + "RF.nc.tif")
    # print(chirps_file)
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

def sifAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,"SIF005_" + str(yy) + str(month).zfill(2) + ".nc.tif")
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file) and os.path.exists(cropland):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -9999) & (crop == 12)
        ave = np.mean(chirps[mask])
        return ave, month
# print(sifAveMap(r"D:\Cornell\EthiopianDrought\AData\SIF005Clip",2015,1))
def sifAveMap2(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,"sif005_eemd_anomaly_" + str(yy) + str(month).zfill(2) + "RF.nc.tif")
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -9999) & (crop == 12)
        ave = np.mean(chirps[mask])
        return ave, month
def calSifAllAve(chirpsdirectory,month):
    # this data from the original data handled by zmh
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "SIF005_" + str(dt.year) + str(dt.month).zfill(2) + ".nc.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(chirpsdirectory,
                                   "SIF005_" + str(dt.year) + str(dt.month).zfill(2) + ".nc.tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = -9999
            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std
def calSifAllAve2(chirpsdirectory,month):
    # this data from the original data handled by zmh
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        chirps_file = os.path.join(chirpsdirectory,
                                   "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = -9999
            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std
def sifAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228, 299))
    multidarr = np.zeros((228, 299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "sif005_eemd_anomaly_" + str(dt.year) + str(dt.month).zfill(2) + "RF.nc.tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = -9999

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1

    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)

    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average, axisTime

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


def eviAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,str(yy)+"."+str(month).zfill(2)+".01"+".tif")
    if not os.path.exists(chirps_file):
        print("can not find the file {}".format(chirps_file))
    # print("*", chirps_file)
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -3000) & (crop ==12)
        ave = np.mean(chirps[mask])
        return ave, month

def calEviAllAve(chirpsdirectory,month):
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
            multidarr[:, :, band_id][crop !=12] = -3000

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -3000)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std

def esacciClip(directory,epregionshppath,ouputdirectory):

    files = glob.glob(os.path.join(directory, "*.tif"))
    # print(chirps_decompress_files)
    for file in files:
        output_raster = os.path.join(ouputdirectory, os.path.basename(file))
        clipbyshp(file, output_raster, epregionshppath, dstNodata=-2)
        print("{} has been processed".format(file))

# esacciClip(r"D:\Cornell\EthiopianDrought\ESACCIverison0.4.5\ESACCITIF",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5")

def calESACCIAverageandStd(directory,outputdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        modis_file = os.path.join(directory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-"+str(dt.year)+str(dt.month).zfill(2) +".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45,59))
    multidarr =np.zeros((45,59,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(directory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")

        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == -2)
            maskarr[mask] =-2
            multidarr[:,:,band_id] = modis
            del modis
            del mask
            band_id += 1


    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis=2)

    data = np.zeros((45, 59, 3))
    data[:, :, 0] = averageMatrix
    data[:, :, 1] = stdMatrix
    data[:, :, 2] = maskarr

    modis_reference =  gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5\ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-199101.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"esacci_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 59, 45, im_bands=3, dtype=gdal.GDT_Float32)
# for i in range(1,13):
#     calESACCIAverageandStd(r"D:\Cornell\EthiopianDrought\AData\ESACCIV0.4.5", r"D:\Cornell\EthiopianDrought\AData\ESACCIPars", i)


def esacciAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45, 59))
    multidarr = np.zeros((45, 59, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(modis == -2)
            # maskarr[mask] = maskarr[mask] + 1

            multidarr[:, :, band_id] = modis
            multidarr[:, :,  band_id][crop != 12] = -2
            axisTime.append(dt.year)
            del modis
            # del mask

            band_id += 1

    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -2)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)

    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average,axisTime

def esacciAnomalyMap(evidirectory,eviparsdirectory,yy,mm):

    evi_file = os.path.join(evidirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(yy) + str(mm).zfill(
                                      2) + ".tif")
    evipars_file = os.path.join(eviparsdirectory,"esacci_month"+str(mm).zfill(2)+".tif")

    evi = gdal.Open(evi_file).ReadAsArray()
    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -2) | (evi == -2)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -2) & (evi > -2)
    mask2 = np.where(mask2)

    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999
    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def esacciAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(yy) + str(month).zfill(
                                      2) + ".tif")
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(yy)))[0]
    # print("*", chirps_file)
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -2) & (crop == 12)
        ave = np.mean(chirps[mask])

        return ave, month


def calEsacciAllAve(chirpsdirectory, month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45, 59))
    multidarr = np.zeros((45, 59, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(modis == -2)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = modis
            multidarr[:,:,band_id][crop != 12] = -2
            axisTime.append(dt.year)
            del modis
            # del mask

            band_id += 1

    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -2)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)


    return mean,month,std


def etClip(directory,epregionshppath,ouputdirectory):

    files = glob.glob(os.path.join(directory, "*.tif"))
    # print(chirps_decompress_files)
    for file in files:
        output_raster = os.path.join(ouputdirectory, os.path.basename(file))
        clipbyshp(file, output_raster, epregionshppath, dstNodata=-999)
        print("{} has been processed".format(file))

# etClip(r"D:\Cornell\EthiopianDrought\AData\ETMonthV3.3a",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\AData\ETClip")

def calETAverageandStd(directory, outputdirectory, month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(directory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45, 59))
    multidarr = np.zeros((45, 59, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(directory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")

        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == -999)
            maskarr[mask] = -999
            multidarr[:, :, band_id] = modis

            del modis
            del mask

            band_id += 1


    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis=2)

    data = np.zeros((45, 59, 3))
    data[:, :, 0] = averageMatrix
    data[:, :, 1] = stdMatrix
    data[:, :, 2] = maskarr
    modis_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\ETClip\ET.v3.3a198101.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory, "et_month" + str(month).zfill(2) + ".tif")
    write_Img(data, outputpath, proj, geotrans, 59, 45, im_bands=3, dtype=gdal.GDT_Float32)


# for i in range(1,13):
#     calETAverageandStd(r"D:\Cornell\EthiopianDrought\AData\ETClip", r"D:\Cornell\EthiopianDrought\AData\ETMonthPars", i)


def etAnomalySeries(chirpsdirectory, month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45, 59))
    multidarr = np.zeros((45, 59, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")

        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(modis == -999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = modis
            multidarr[:,:,band_id][crop != 12] = -999
            axisTime.append(dt.year)
            del modis
            # del mask

            band_id += 1

    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average, axisTime


def etAnomalyMap(evidirectory, eviparsdirectory, yy, mm):
    evi_file = os.path.join(evidirectory,
                            "ET.v3.3a" + str(yy) + str(mm).zfill(
                                2) + ".tif")
    evipars_file = os.path.join(eviparsdirectory, "et_month" + str(mm).zfill(2) + ".tif")

    evi = gdal.Open(evi_file).ReadAsArray()

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -999) | (evi == -999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -999) & (evi > -999)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999

    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def etAveMap(chirpsdirectory,yy, month):


    chirps_file =os.path.join(chirpsdirectory,
                            "ET.v3.3a" + str(yy) + str(month).zfill(
                                2) + ".tif")
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(yy)))[0]
    # print("*", chirps_file)
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -999) & (crop == 12)
        ave = np.mean(chirps[mask])

        return ave, month

def calEtAllAve(chirpsdirectory, month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((45, 59))
    multidarr = np.zeros((45, 59, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(chirpsdirectory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1v006AggClip', "*MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(modis == -999)
            # maskarr[mask] = maskarr[mask] + 1
            multidarr[:, :, band_id] = modis
            multidarr[:,:,band_id][crop != 12] =-999
            axisTime.append(dt.year)
            del modis
            # del mask

            band_id += 1


    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std

def modisETclip(modisdirectory,epregionshppath,clippeddirectory):

    modis_files = glob.glob(os.path.join(modisdirectory,"*"))
    # print(chirps_decompress_files)
    for subfile in modis_files:

        input_raster = subfile
        output_raster = os.path.join(clippeddirectory,os.path.basename(subfile)+".tif")


        # gdal.Warp(r"D:\Cornell\EthiopianDrought\Test\2000.02.01all.tif", subdateset, dstSRS='EPSG:4326',
        #           dstNodata=-9999)

        clipbyshp(input_raster, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=255)
        print("{} has been processed".format(input_raster))

# modisETclip(r"D:\Cornell\BessMonth",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\BessMonthClip")


def calModisETAverageandStd(modisdirectory,outputdirectory,month):
    start = datetime.strptime("-".join(["2003",str(month).zfill(2),"01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,
                                  str(dt.year) + str(dt.month).zfill(2) + str(dt.day).zfill(2) + ".tif.tif")
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((1378,1779))
    multidarr =np.zeros((1378,1779,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,str(dt.year) +str(dt.month).zfill(2)+str(dt.day).zfill(2)+".tif.tif")

        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == 255)
            maskarr[mask] = 255
            multidarr[:,:,band_id] = modis
            del modis
            del mask
            band_id += 1

    averageMatrix = multidarr.mean(axis=2)
    stdMatrix = multidarr.std(axis= 2)

    data = np.zeros((1378,1779,3))
    data[:,:,0] = averageMatrix
    data[:,:,1] = stdMatrix
    data[:,:,2] = maskarr

    modis_reference =  gdal.Open(r"D:\Cornell\BessMonthClip\20000301.tif.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"bess_et_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 1779, 1378, im_bands=3, dtype=gdal.GDT_Float32)
# for i in range(1,13):
#     calModisETAverageandStd(r"D:\Cornell\BessMonthClip",
#                      r"D:\Cornell\EthiopianDrought\AData\BessETMonthPars", i)

def BessETAnomalyMap(evidirectory,eviparsdirectory,yy,mm):

    evi_file = os.path.join(evidirectory,str(yy)+str(mm).zfill(2)+"01.tif.tif")
    evipars_file = os.path.join(eviparsdirectory,"bess_et_month"+str(mm).zfill(2)+".tif")

    evi = gdal.Open(evi_file).ReadAsArray()

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == 255) | (evi == 255)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr < 255) & (evi < 255)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999
    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def BessETAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year)  + str(dt.month).zfill(2) +"01" + ".tif.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((1378,1779))
    multidarr = np.zeros((1378,1779, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year) + str(dt.month).zfill(2) +"01" + ".tif.tif")


        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12Q1AggClip', "agg_{}*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = 255

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband < 255)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average,axisTime


def BessETAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,str(yy)+str(month).zfill(2)+"01"+".tif.tif")
    if not os.path.exists(chirps_file):
        print("can not find the file {}".format(chirps_file))
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12Q1AggClip', "agg_{}*.tif".format(dt.year)))[0]
    # print("*", chirps_file)
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps < 255) & (crop ==12)
        ave = np.mean(chirps[mask])
        return ave, month

def calBessETAllAve(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year) + str(dt.month).zfill(2) +"01" + ".tif.tif")
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((1378,1779))
    multidarr = np.zeros((1378,1779, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   str(dt.year) + str(dt.month).zfill(2) +"01" + ".tif.tif")
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12Q1AggClip', "agg_{}*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:, :, band_id][crop != 12] = 255
            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband < 255)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std

def GoSIFclip(modisdirectory,epregionshppath,clippeddirectory):

    modis_files = glob.glob(os.path.join(modisdirectory,"*"))
    # print(chirps_decompress_files)
    for subfile in modis_files:

        input_raster = subfile
        print(input_raster)
        output_raster = os.path.join(clippeddirectory,os.path.basename(subfile))


        # gdal.Warp(r"D:\Cornell\EthiopianDrought\Test\2000.02.01all.tif", subdateset, dstSRS='EPSG:4326',
        #           dstNodata=-9999)

        clipbyshp(input_raster, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=32767)
        print("{} has been processed".format(input_raster))

# GoSIFclip(r"D:\Cornell\GosifMonthly",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\GOSIFClip")


# GoSIFclip(r"D:\Cornell\GOSIFV002",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\GOSIFV002Clip")

def calGOSIFAverageandStd(modisdirectory,outputdirectory,month):
    start = datetime.strptime("-".join(["2003",str(month).zfill(2),"01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,"{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr =np.zeros((228,299,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,"{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))

        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == -9999)
            maskarr[mask] = -9999
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

    modis_reference =  gdal.Open(r"D:\Cornell\GOSIFV002Clip\GOSIF_2000.M03.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"GOSIF_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)
# for i in range(1,13):
#     calGOSIFAverageandStd(r"D:\Cornell\EthiopianDrought\AData\GOSIF0317AnomSplit",
#                      r"D:\Cornell\EthiopianDrought\AData\GOSIFMonthPars", i)

def GOSIFAnomalyMap(evidirectory,eviparsdirectory,yy,mm):

    evi_file = os.path.join(evidirectory,"{}GOSIF_Anom{}.tif".format(str(yy),str(mm).zfill(2)))
    evipars_file = os.path.join(eviparsdirectory,"GOSIF_month"+str(mm).zfill(2)+".tif")
    print(evi_file)
    evi = gdal.Open(evi_file).ReadAsArray()

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (evi == -9999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (evi > -9999)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999
    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def GOSIFAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr = np.zeros((228,299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:,:,band_id][crop !=12] = -9999

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average,axisTime


def GOSIFAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,"{}GOSIF_Anom{}.tif".format(str(yy),str(month).zfill(2)))
    if not os.path.exists(chirps_file):
        print("can not find the file {}".format(chirps_file))
    # print("*", chirps_file)
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -9999) & (crop == 12)
        ave = np.mean(chirps[mask])
        return ave, month

def calGOSIFAllAve(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr = np.zeros((228,299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}GOSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        # print("*", chirps_file)
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:,:,band_id][crop != 12] = -9999

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std
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




def calNewSIFAverageandStd(modisdirectory,outputdirectory,month):
    start = datetime.strptime("-".join(["2007",str(month).zfill(2),"01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(modis_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr =np.zeros((228,299,bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        modis_file = os.path.join(modisdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))

        if os.path.exists(modis_file):
            modis = gdal.Open(modis_file).ReadAsArray()
            mask = np.where(modis == -9999)
            maskarr[mask] = -9999
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

    modis_reference =  gdal.Open(r"D:\Cornell\NewSIF005Clip\SIF005_201512.nc.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"NewSIF_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)
# for i in range(1,13):
#     calNewSIFAverageandStd(r"D:\Cornell\EthiopianDrought\AData\NewSIF0318AnomSplit2007",
#                      r"D:\Cornell\EthiopianDrought\AData\NewSIFMonthPars2007", i)

def NewSIFAnomalyMap(evidirectory,eviparsdirectory,yy,mm):

    evi_file = os.path.join(evidirectory,"{}NewSIF_Anom{}.tif".format(str(yy),str(mm).zfill(2)))
    evipars_file = os.path.join(eviparsdirectory,"NewSIF_month"+str(mm).zfill(2)+".tif")

    evi = gdal.Open(evi_file).ReadAsArray()

    evipars = gdal.Open(evipars_file)
    averageMatrix = evipars.GetRasterBand(1).ReadAsArray()
    stdMatrix = evipars.GetRasterBand(2).ReadAsArray()
    maskarr = evipars.GetRasterBand(3).ReadAsArray()

    mask = (stdMatrix <= 0) | (maskarr == -9999) | (evi == -9999)
    mask = np.where(mask)
    mask2 = (stdMatrix > 0) & (maskarr > -9999) & (evi > -9999)
    mask2 = np.where(mask2)
    anomalyMap = np.zeros(shape=evi.shape)
    anomalyMap[mask] = -9999
    anomalyMap[mask2] = (evi[mask2] - averageMatrix[mask2])/stdMatrix[mask2]

    return anomalyMap

def NewSIFAnomalySeries(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2007", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,"{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr = np.zeros((228,299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        # print("*", chirps_file)
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:,:,band_id][crop !=12] = -9999

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []

    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)
    average = (ave - mean)/std
    return average,axisTime


def NewSIFAveMap(chirpsdirectory,yy, month):

    chirps_file =os.path.join(chirpsdirectory,"{}NewSIF_Anom{}.tif".format(str(yy),str(month).zfill(2)))
    if not os.path.exists(chirps_file):
        print("can not find the file {}".format(chirps_file))
    # print("*", chirps_file)
    cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(yy)))[0]
    if os.path.exists(chirps_file):
        chirps = gdal.Open(chirps_file).ReadAsArray()
        crop = gdal.Open(cropland).ReadAsArray()
        mask = (chirps > -9999) & (crop == 12)
        ave = np.mean(chirps[mask])
        return ave, month

def calNewSIFAllAve(chirpsdirectory,month):
    start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    bandNum = 0
    axisTime = []
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        if os.path.exists(chirps_file):
            bandNum += 1

    maskarr = np.zeros((228,299))
    multidarr = np.zeros((228,299, bandNum))
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        chirps_file = os.path.join(chirpsdirectory,
                                   "{}NewSIF_Anom{}.tif".format(str(dt.year),str(dt.month).zfill(2)))
        # print("*", chirps_file)
        cropland = glob.glob(os.path.join(r'D:\Cornell\MCD12C1V006Clip', "MCD12C1.A{}001.*.tif".format(dt.year)))[0]
        if os.path.exists(chirps_file):
            chirps = gdal.Open(chirps_file).ReadAsArray()
            crop = gdal.Open(cropland).ReadAsArray()
            # mask = np.where(chirps == -9999)
            # maskarr[mask] = -9999
            multidarr[:, :, band_id] = chirps
            multidarr[:,:,band_id][crop != 12] = -9999

            del chirps
            del crop
            # del mask
            axisTime.append(dt.year)
            band_id += 1


    ave = []
    num = []
    for band in range(bandNum):
        mband = multidarr[:, :, band]
        mask = np.where(mband > -9999)
        ave.append(np.mean(mband[mask]))

    ave = np.array(ave)
    mean = np.mean(ave)
    std = np.std(ave)

    return mean,month,std