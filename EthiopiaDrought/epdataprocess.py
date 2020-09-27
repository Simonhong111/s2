from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
def clipbyshp(input_raster,output_raster,input_shape):
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
                   dstNodata=-9999)  # select the no data value you like
    ds = None

def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)

    dataset.SetProjection(str(proj))
    for id in range(im_bands):
        print("**********")
        dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset

def chirpsclip(chirpsdirectory,epregionshppath,clippeddirectory):

    chirps_decompress_files = glob.glob(os.path.join(chirpsdirectory,"*.tif"))
    # print(chirps_decompress_files)
    for tiffile in chirps_decompress_files:
        input_raster = glob.glob(os.path.join(tiffile,"*.tif"))[0]
        output_raster = os.path.join(clippeddirectory,os.path.basename(input_raster))
        clipbyshp(input_raster,output_raster,epregionshppath)
        print("{} has been processed".format(input_raster))


# chirpsclip(r"D:\Cornell\EthiopianDrought\tifs",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               "D:\Cornell\EthiopianDrought\Chirps2")

def calChirpsAverageandStd(chirpsclippeddirectory,outputdirectory,month):


        start = datetime.strptime("-".join(["1981", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
        stop = datetime.strptime("2020-01-01", "%Y-%m-%d").date()
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



def chirpsAnomaly(chirpsaveragepath,mm):

    chirpsaverage = pd.read_csv(chirpsaveragepath)
    # print(chirpsaverage[0:10])
    ave_all = chirpsaverage.average[chirpsaverage.month==mm].mean()
    std_all = chirpsaverage.average[chirpsaverage.month==mm].std()
    year = chirpsaverage.year[chirpsaverage.month==mm]
    anomaly = (chirpsaverage.average[chirpsaverage.month==mm]- ave_all)/std_all
    # print(len(year),ave_all,std_all,anomaly)
    return anomaly,year
# anomaly,year  = chirpsAnomaly(r"D:\Cornell\EthiopianDrought\Chirps2CSV\EthiopiaRegion.csv",1)


def generateSIFAnomaly(sifanomalypath,sifouputdirectory):

    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Test\chirps-v2.0.1981.01.tif")
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
    mask = np.where(sifanomaly < -9999)
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
def sifaverage(sifclippeddirectory,sifaverageoutputpath):

    sifclippedfiles = glob.glob(os.path.join(sifclippeddirectory,"*.tif"))
    Ave = []
    YY = []
    Month = []
    for clippedfile in sifclippedfiles:

        print(clippedfile)
        basename = os.path.basename(clippedfile)
        year,month = basename.split("_")[3][0:4],basename.split("_")[3][4:6]
        raster = gdal.Open(clippedfile,0).ReadAsArray()
        mask = np.where(raster >-9999.0)
        average = np.mean(raster[mask])
        Ave.append(average)
        YY.append(year)
        Month.append(month)
    chirpsdict = {"year":YY,"month":Month,"average":Ave}
    chirpsPd = pd.DataFrame(chirpsdict)
    chirpsPd.to_csv(sifaverageoutputpath)
# sifaverage(r"D:\Cornell\EthiopianDrought\clipedSIf",r"D:\Cornell\EthiopianDrought\SIFCSV\EthiopiaRegionSIFAverage.csv")
def sifAnomaly2(sifaveragepath,mm):

    sifaverage = pd.read_csv(sifaveragepath)
    # print(chirpsaverage[0:10])
    ave_all = sifaverage.average[sifaverage.month==mm]
    year = sifaverage.year[sifaverage.month==mm]
    print("*",ave_all)
    # print(len(year),ave_all,std_all,anomaly)
    return ave_all,year

# def generateSIFAnomaly2(sifanomalypath,sifouputdirectory):
#
#     raster = gdal.Open(r"D:\Cornell\EthiopianDrought\Test\chirps-v2.0.1981.01.tif")
#     orgLon = -180
#     orgLat = 90
#
#     geotrans = [orgLon,0.05,0.0,orgLat,0.0,-0.05]
#     proj = raster.GetProjection()
#     sifoutputpath = os.path.join(sifouputdirectory,os.path.basename(sifanomalypath)+"all.tif")
#
#     fin = Dataset(sifanomalypath,"r")
#     sifanomaly = fin.variables['SIF_740_daily_corr_anomaly'][:]
#
#     write_Img(sifanomaly, sifoutputpath, proj, geotrans, 7200, 3600, im_bands=1, dtype=gdal.GDT_Float32)
#     fin.close()
#     del sifanomaly
#
# generateSIFAnomaly2(r"D:\Cornell\EthiopianDrought\sif005_eemd_anomaly_200211RF.nc",r"D:\Cornell\EthiopianDrought\Test")





# generateSIFAnomaly(r"D:\Cornell\EthiopianDrought\sif005_eemd_anomaly_200211RF.nc",r"D:\Cornell\EthiopianDrought\Test")

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

        clipbyshp(subdateset, output_raster, r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp")
        print("{} has been processed".format(input_raster))



# modisclip(r"D:\Cornell\EthiopianDrought\MOD13C2.006",
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
#               r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia")


def calModisAverageandStd(modisdirectory,outputdirectory,month):

    start = datetime.strptime("-".join(["2000",str(month).zfill(2),"01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2020-01-01", "%Y-%m-%d").date()
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
        print("*",modis_file)
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

    modis_reference =  gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    geotrans = modis_reference.GetGeoTransform()
    proj = modis_reference.GetProjection()
    outputpath = os.path.join(outputdirectory,"evi_month"+str(month).zfill(2)+".tif")
    write_Img(data,outputpath , proj, geotrans, 299, 228, im_bands=3, dtype=gdal.GDT_Float32)



#
# for i in range(1,13):
#     calModisAverageandStd(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia",
#                      "D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyPars", i)

# modis_files = glob.glob(os.path.join(r"D:\Cornell\EthiopianDrought\Test\modis","*"))
#     # print(chirps_decompress_files)
# for subfile in modis_files:
#
#     input_raster = glob.glob(os.path.join(subfile,"*.hdf"))[0]
#     output_raster = os.path.join(r"D:\Cornell\EthiopianDrought\Test",os.path.basename(subfile)+".tif")
#
#     modis = gdal.Open(input_raster)
#     subdateset = modis.GetSubDatasets()[1][0]
#     print(subdateset)
#     gdal.Warp(r"D:\Cornell\EthiopianDrought\Test\2000.02.01all.tif", subdateset, dstSRS='EPSG:4326',dstNodata=-9999)
#
#     clipbyshp(subdateset,output_raster,r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp")
#     print("{} has been processed".format(input_raster))


