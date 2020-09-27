from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
from netCDF4 import Dataset
from Experiment.Exp_Clip_Raster_By_SHP import clipbyshp
from Experiment.Exp_Write_Image import write_Img
def GenerateNewSIF(SIFPath,OutDir):

    raster = gdal.Open(r"D:\Cornell\EthiopianDrought\SIFAnomaly\sif005_eemd_anomaly_200208RF.nc.tif")
    orgLon = 32.8
    orgLat = 15
    RasterXSize = 310
    RasterYSize = 240
    geotrans = [orgLon,0.05,0.0,orgLat,0.0,-0.05]
    proj = raster.GetProjection()

    sifoutputpath = os.path.join(OutDir,os.path.basename(SIFPath)+".tif")
    eLon = orgLon + RasterXSize*0.05
    eLat = orgLat - RasterYSize*0.05

    orgLon_Col = int((orgLon + 180)*20)
    orgLat_Row = int((90 - orgLat)*20)

    eLon_Col = int((eLon + 180)*20)
    eLat_Row = int((90 - eLat)*20)

    fin = Dataset(SIFPath,"r")
    sif_data = fin.variables['SIF_740_daily_corr'][:]
    mask = np.where(sif_data <= -999)
    sif_data[mask] = -9999
    sif = sif_data[orgLat_Row:eLat_Row,orgLon_Col:eLon_Col]
    write_Img(sif, sifoutputpath, proj, geotrans, 310, 240, im_bands=1, dtype=gdal.GDT_Float32)
    fin.close()
    del sif_data
    del sif

# sif__files = glob.glob(os.path.join(r"D:\Cornell\JiaMing\SIFV2\SIF005v20200903", "*.nc"))
# for sif_path in sif__files:
#     GenerateNewSIF(sif_path, r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_OutRect_Ethiopia")

def SIFClip(SIFDir,epregionshppath,sifouputdirectory):

    sif__files = glob.glob(os.path.join(SIFDir, "*.tif"))

    for sif_path in sif__files:
        output_raster = os.path.join(sifouputdirectory, os.path.basename(sif_path))
        clipbyshp(sif_path, output_raster, epregionshppath,dstNodata=-9999)
        print("{} has been processed".format(sif_path))

SIFClip(r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_OutRect_Ethiopia",
              r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp",
              r"D:\Cornell\EthiopianDrought\0ExperimentData\SIF_Data\SIF_V3_Ethiopia")

# path = r"D:\Cornell\JiaMing\SIFV2\SIF005v20200903\SIF005_200208.nc"
# out =r"D:\Cornell\EthiopianDrought\0ExperimentData\TempData\SIF005_200208.nc.tif"
# fin = Dataset(path,"r")
# sif_data = fin.variables['SIF_740_daily_corr'][:]
# mask = np.where(sif_data <= -999)
# sif_data[mask] = -9999
# geotrans = [-180,0.05,0,90,0,-0.05]
# proj = osr.SpatialReference()
# proj.ImportFromEPSG(4326)
# write_Img(sif_data, out, proj, geotrans, 7200, 3600, im_bands=1, dtype=gdal.GDT_Float32)
# fin.close()
