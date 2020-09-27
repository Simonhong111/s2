from osgeo import gdal,osr,ogr
import os
import numpy as np
from dateutil import rrule
from datetime import *
from scipy import signal
from Experiment.Exp_Write_Image import write_Img

def fit(y):
    return signal.detrend(y)

def detrendSIF(SIFDir,OutDir,month):
    st_year = 2007
    start = datetime.strptime("-".join([str(st_year), str(month).zfill(2), "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()

    refpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia\SIF005_200208.nc.tif"
    reference = gdal.Open(refpath)
    geotrans = reference.GetGeoTransform()
    proj = reference.GetProjection()
    Width, Height = reference.RasterXSize, reference.RasterYSize

    bandNum = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        sif_file = os.path.join(SIFDir,
                                   "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))
        if os.path.exists(sif_file):
            bandNum += 1
        assert os.path.exists(sif_file),"{} does not exist".format(sif_file)

    maskarr = np.zeros((Height, Width))
    multidarr = np.zeros((Height, Width, bandNum),dtype=np.float)
    band_id = 0

    for dt in (rrule.rrule(rrule.YEARLY, interval=1, dtstart=start, until=stop)):

        sif_file = os.path.join(SIFDir,
                                   "SIF005_{}{}.nc.tif".format(str(dt.year), str(dt.month).zfill(2)))

        if os.path.exists(sif_file):
            sif_data = gdal.Open(sif_file).ReadAsArray()
            maskarr[sif_data < -0.05] = -9999
            multidarr[:, :, band_id] = sif_data
            band_id += 1

    multiarrflatten = multidarr.reshape(Width*Height,bandNum)
    results = map(fit,multiarrflatten)
    deResults = np.array(list(results)).reshape(Height,Width,bandNum)

    for i in range(bandNum):
        deResults[:,:,i][maskarr ==-9999] = -9999
        outputpath = os.path.join(OutDir, "Detrend_SIF005_{}{}.nc.tif".format(i +st_year,str(month).zfill(2)))
        write_Img(deResults[:,:,i], outputpath, proj, geotrans, Width, Height, im_bands=1, dtype=gdal.GDT_Float32)


for i in range(1,13):
    res = detrendSIF(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia",
                     r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\Detrend_SIF_V3_Ethiopia",i)




