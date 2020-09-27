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
from matplotlib import pyplot as plt
import csv
from scipy import stats
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

# path = r"D:\Cornell\EthiopianDrought\CropType2015\GLASS-GLC_7classes_2015.tif"
# output = r"D:\Cornell\EthiopianDrought\CropType2015\GLASS-landcover_2015.tif"
# clipbyshp(path,
#               output,
#               r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",dstNodata=0)


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

def write2csv(PVIlongDir,PVIshortDir,EVIDir,RainFallDir):

    landcover = gdal.Open(r"D:\Cornell\EthiopianDrought\CropType2015\GLASS-landcover_2015.tif")
    ldproj = landcover.GetProjection()
    ldgeo = landcover.GetGeoTransform()
    ldproj_new = osr.SpatialReference()
    ldproj_new.ImportFromWkt(str(ldproj))

    PVI = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\PVI\long_pvi_2010.tif")
    pproj = PVI.GetProjection()
    pgeo = PVI.GetGeoTransform()
    pproj_new = osr.SpatialReference()
    pproj_new.ImportFromWkt(str(pproj))
    pWidth = PVI.RasterXSize

    EVI = gdal.Open(r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia\2000.02.01.tif")
    eproj = EVI.GetProjection()
    egeo = EVI.GetGeoTransform()
    eproj_new = osr.SpatialReference()
    eproj_new.ImportFromWkt(str(eproj))
    eWidth = EVI.RasterXSize

    RF = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.1981.01.tif")
    rproj = RF.GetProjection()
    rgeo = RF.GetGeoTransform()
    rproj_new = osr.SpatialReference()
    rproj_new.ImportFromWkt(str(rproj))
    rWidth = RF.RasterXSize

    ld_pvi_ct = osr.CoordinateTransformation(ldproj_new,pproj_new)
    ld_evi_ct = osr.CoordinateTransformation(ldproj_new,eproj_new)
    ld_rf_ct = osr.CoordinateTransformation(ldproj_new,rproj_new)

    landraster = landcover.ReadAsArray()
    mask = landraster == 10
    cropROW,cropCOL = np.where(landraster == 10)
    cropCOL = cropCOL + 0.5
    cropROW = cropROW + 0.5
    print("hh",list(zip(cropROW,cropCOL)))
    cropEast = ldgeo[0] + cropCOL*ldgeo[1]
    cropNorth = ldgeo[3] + cropROW*ldgeo[5]

    cropPoints = np.zeros(shape=(len(cropEast),2),dtype=np.float64)
    cropPoints[:,0] = cropEast
    cropPoints[:,1] = cropNorth

    pviPoints = ld_pvi_ct.TransformPoints(cropPoints)
    pviPoints = np.array(pviPoints)
    pviCOL = np.array((pviPoints[:,0] - pgeo[0])/pgeo[1]).astype(np.int)
    pviROW = np.array((pviPoints[:, 1] - pgeo[3]) / pgeo[5]).astype(np.int)
    pInd = pviROW.flatten()*pWidth + pviCOL.flatten()


    eviPoints = ld_evi_ct.TransformPoints(cropPoints)
    eviPoints = np.array(eviPoints)
    eviCOL = np.array((eviPoints[:, 0] - egeo[0]) / egeo[1]).astype(np.int)
    eviROW = np.array((eviPoints[:, 1] - egeo[3]) / egeo[5]).astype(np.int)
    eInd = eviROW.flatten() * eWidth + eviCOL.flatten()



    rfPoints = ld_rf_ct.TransformPoints(cropPoints)
    rfPoints = np.array(rfPoints)
    rfCOL = np.array((rfPoints[:, 0] - rgeo[0]) / rgeo[1]).astype(np.int)
    rfROW = np.array((rfPoints[:, 1] - rgeo[3]) / rgeo[5]).astype(np.int)
    rInd = rfROW.flatten() * rWidth + rfCOL.flatten()

    PVIshort = []
    PVIlong = []
    EVIData = []
    RFData = []

    for month in range(2,10):

        start = datetime.strptime("-".join(["2003", str(month).zfill(2), "01"]), "%Y-%m-%d").date()
        stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
        for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
            evi_file = os.path.join(EVIDir,"{}.{}.01.tif".format(str(dt.year), str(dt.month).zfill(2)))
            assert os.path.exists(evi_file),"the file {} does not exist".format(evi_file)

            rf_file = os.path.join(RainFallDir, "chirps-v2.0.{}.{}.tif".format(str(dt.year), str(dt.month).zfill(2)))
            assert os.path.exists(rf_file), "the file {} does not exist".format(rf_file)

            eviraster = gdal.Open(evi_file).ReadAsArray()
            evicrop = np.take(eviraster,eInd)
            EVIData.append(evicrop.tolist())

            rfraster = gdal.Open(rf_file).ReadAsArray()
            rfcrop = np.take(rfraster,rInd)
            RFData.append(rfcrop.tolist())

    start = datetime.strptime("2003-01-01", "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):
        long_file = os.path.join(PVIlongDir, "long_pvi_{}.tif".format(str(dt.year)))
        assert os.path.exists(long_file), "the file {} does not exist".format(long_file)

        short_file = os.path.join(PVIshortDir, "short_pvi_{}.tif".format(str(dt.year)))
        assert os.path.exists(short_file), "the file {} does not exist".format(short_file)

        shortraster = gdal.Open(short_file).ReadAsArray()
        shortcrop = np.take(shortraster, pInd)
        PVIshort.append(shortcrop.tolist())

        longraster = gdal.Open(long_file).ReadAsArray()
        longcrop = np.take(longraster, pInd)
        PVIlong.append(longcrop.tolist())

    PVIshort = np.array(PVIshort)
    PVIlong = np.array(PVIlong)
    EVIData = np.array(EVIData)
    RFData = np.array(RFData)

    DataDict = {}
    print(PVIlong.shape,PVIshort.shape,EVIData.shape,RFData.shape)
    index = 0
    for month in range(2,10):
        for year in range(2003,2019):

            DataDict["EVI"+str(year)+str(month).zfill(2)] = EVIData[index,:]

            index += 1
    index = 0
    for month in range(2,10):
        for year in range(2003,2019):


            DataDict["RF"+str(year)+str(month).zfill(2)] = RFData[index,:]
            index += 1
    index2 = 0
    for year in range(2003,2019):
        DataDict["PVIshort" + str(year)] = PVIshort[index2, :]

        index2 += 1
    index2 = 0
    for year in range(2003, 2019):

        DataDict["PVIlong" + str(year)] = PVIlong[index2, :]
        index2 += 1
    # df = pd.DataFrame(DataDict)
    # outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Variable_2003-2018.csv"
    # df.to_csv(outpath, index=False)





    # plt.imshow(landraster)
    # plt.imshow(mask*1)
    # plt.show()


PVIlongDir = r"D:\Cornell\EthiopianDrought\AData\PVI"
PVIshortDir =r"D:\Cornell\EthiopianDrought\AData\PVI"
EVIDir = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
RainFallDir = r"D:\Cornell\EthiopianDrought\Chirps2"

# write2csv(PVIlongDir,PVIshortDir,EVIDir,RainFallDir)

def valid():
    varpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Variable_2003-2018.csv"
    Titles = []
    for month in range(2, 10):
        for year in range(2003, 2019):
            Titles.append("EVI" + str(year) + str(month).zfill(2))

    for month in range(2, 10):
        for year in range(2003, 2019):
            Titles.append("RF" + str(year) + str(month).zfill(2))

    for year in range(2003, 2019):
        Titles.append("PVIshort" + str(year))

    for year in range(2003, 2019):
        Titles.append("PVIlong" + str(year))

    VData = []
    with open(varpath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:

            temp = []
            areatemp = []
            for i in range(288):
                temp.append(float(row[Titles[i]]))

            if np.array(temp).min() <= -3000:
                continue
            VData.append(temp)
    VData = np.array(VData)
    print(VData.shape, "sh")
    VDict = {}
    for i in range(288):
        VDict[Titles[i]] = VData[:, i]

    # df = pd.DataFrame(VDict)
    # outpath = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Variable_2003-2018_Valid.csv"
    # df.to_csv(outpath, index=False)
    return VData




VData =  valid()
def extract(VData):
    path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\Variable_2003-2018_Valid.csv"
    SEVI = []
    LEVI = []
    PVIShort =[]
    PVILong = []
    SRF =[]
    LRF = []

    for year in range(2003,2019):
        shortevisub = []
        longevitsub = []
        shortrftsub = []
        longrftsub = []
        shortpvisub =[]
        longpvitsub =[]
        with open(path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                shortevitemp = []
                longevitemp =[]
                shortrftemp =[]
                longrftemp=[]
                shortpvitemp = []
                longpvitemp =[]
                shortevitemp.append(float(row["EVI"+str(year)+"02"]))
                shortevitemp.append(float(row["EVI" + str(year) + "03"]))
                shortevitemp.append(float(row["EVI" + str(year) + "04"]))
                shortevitemp.append(float(row["EVI" + str(year) + "05"]))
                longevitemp.append(float(row["EVI" + str(year) + "06"]))
                longevitemp.append(float(row["EVI" + str(year) + "07"]))
                longevitemp.append(float(row["EVI" + str(year) + "08"]))
                longevitemp.append(float(row["EVI" + str(year) + "09"]))
                shortrftemp.append(float(row["RF" + str(year) + "02"]))
                shortrftemp.append(float(row["RF" + str(year) + "03"]))
                shortrftemp.append(float(row["RF" + str(year) + "04"]))
                shortrftemp.append(float(row["RF" + str(year) + "05"]))
                longrftemp.append(float(row["RF" + str(year) + "06"]))
                longrftemp.append(float(row["RF" + str(year) + "07"]))
                longrftemp.append(float(row["RF" + str(year) + "08"]))
                longrftemp.append(float(row["RF" + str(year) + "09"]))

                shortpvitemp.append(float(row["PVIshort" + str(year)]))
                longpvitemp.append(float(row["PVIlong" + str(year)]))

                shortevisub.append(np.array(shortevitemp).mean())
                longevitsub.append(np.array(longevitemp).mean())
                shortrftsub.append(np.array(shortrftemp).mean())
                longrftsub.append(np.array(longrftemp).mean())
                shortpvisub.append(np.array(shortpvitemp).mean())
                longpvitsub.append(np.array(longpvitemp).mean())


        SEVI.append(shortevisub)
        LEVI.append(longevitsub)
        SRF.append(shortrftsub)
        LRF.append(longrftsub)
        PVIShort.append(shortpvisub)
        PVILong.append(longpvitsub)

    SEVI = np.array(SEVI)
    LEVI = np.array(LEVI)
    SRF = np.array(SRF)
    LRF = np.array(LRF)
    PVIShort = np.array(PVIShort)
    PVILong = np.array(PVILong)
    print(SEVI.shape,LEVI.shape,SRF.shape,LRF.shape,PVIShort.shape,PVILong.shape)

    SEVIMean = SEVI.mean(axis=0)
    SEVIStd = SEVI.std(axis=0)
    SEVIAnom = np.zeros_like(SEVI,dtype=np.float64)
    for i in range(450):
        SEVIAnom[:,i] = (SEVI[:,i] - SEVIMean[i] )/SEVIStd[i]

    LEVIMean = LEVI.mean(axis=0)
    LEVIStd = LEVI.std(axis=0)
    LEVIAnom = np.zeros_like(LEVI, dtype=np.float64)
    for i in range(450):
        LEVIAnom[:, i] = (LEVI[:, i] - LEVIMean[i]) / LEVIStd[i]

    SRFMean = SRF.mean(axis=0)
    SRFStd = SRF.std(axis=0)
    SRFAnom = np.zeros_like(SRF, dtype=np.float64)
    for i in range(450):
        SRFAnom[:, i] = (SRF[:, i] - SRFMean[i]) / SRFStd[i]
    print(SRFMean.shape,"**")
    LRFMean = LRF.mean(axis=0)
    LRFStd = LRF.std(axis=0)
    LRFAnom = np.zeros_like(LRF, dtype=np.float64)
    for i in range(450):
        LRFAnom[:, i] = (LRF[:, i] - LRFMean[i]) / LRFStd[i]

    SPVIMean = PVIShort.mean(axis=0)
    SPVIStd = PVIShort.std(axis=0)
    SPVIAnom = np.zeros_like(PVIShort, dtype=np.float64)
    for i in range(450):
        SPVIAnom[:, i] = (PVIShort[:, i] - SPVIMean[i]) / SPVIStd[i]

    LPVIMean = PVILong.mean(axis=0)
    LPVIStd = PVILong.std(axis=0)
    LPVIAnom = np.zeros_like(PVILong, dtype=np.float64)
    for i in range(450):
        LPVIAnom[:, i] = (PVILong[:, i] - LPVIMean[i]) / LPVIStd[i]

    SCorrPVIEVI = []
    LCorrPVIEVI = []
    SCorrRFEVI = []
    LCorrRFEVI = []
    for year in range(16):
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(SEVI[year,:],
                                                                       PVIShort[year,:])
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(LEVI[year,:],
                                                                            PVILong[year,:])

        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(SEVI[year, :],
                                                                            SRF[year, :])
        slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(LEVI[year, :],
                                                                            LRF[year, :])
        SCorrPVIEVI.append(r_value1*(-1))
        LCorrPVIEVI.append(r_value2*(-1))
        SCorrRFEVI.append(r_value3)
        LCorrRFEVI.append(r_value4)
    return SCorrPVIEVI,LCorrPVIEVI,SCorrRFEVI,LCorrRFEVI,



SCorrPVIEVI,LCorrPVIEVI,SCorrRFEVI,LCorrRFEVI = extract(VData)
year = range(2003,2019)
plt.xlabel("year")
plt.ylabel("Correlation R")
plt.plot(year,SCorrPVIEVI,label="original short rain pvi evi")
# plt.plot(year,LCorrPVIEVI,label="original long rain pvi evi")
plt.plot(year,SCorrRFEVI,label="original short rain rain evi")
# plt.plot(year,LCorrRFEVI,label="original long rain rain evi")
plt.legend()
plt.show()


