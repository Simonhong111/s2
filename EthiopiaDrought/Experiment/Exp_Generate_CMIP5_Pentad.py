from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
from calendar import monthrange
from Experiment.Exp_Write_Image import write_Img as Write
def GenCmip5Pentand(cmip5Dir,OutDir,yy,mm,pentandID):
    RefPath = r"D:\Cornell\EthiopianDrought\CMIP5DailyClip\Big\cmip5_20060101.tif"
    RefImg = gdal.Open(RefPath)
    RefGeoTrans = RefImg.GetGeoTransform()
    RefProj = RefImg.GetProjection()
    RefW, RefH = RefImg.RasterXSize, RefImg.RasterYSize

    if pentandID <=5:
        DayNum = 5
    if pentandID ==6:
        DayNum = monthrange(yy,mm)[1] - 25

    DaysofMonth =  monthrange(yy,mm)[1]

    if pentandID <=5:
        st_day = (pentandID-1)*5+1
        en_day = st_day+5
    if pentandID == 6:
        st_day = (pentandID - 1) * 5 + 1
        en_day = DaysofMonth + 1

    PenIMGs = np.zeros(shape=(RefH, RefW, DayNum), dtype=np.float)
    PenMask = np.zeros(shape=(RefH, RefW))
    for id,dayId in enumerate(range(st_day,en_day)):
        cmip_file = os.path.join(cmip5Dir,"cmip5_{}{}{}.tif").format(str(yy),str(mm).zfill(2),str(dayId).zfill(2))
        print(cmip_file,dayId)
        cmip_data = gdal.Open(cmip_file).ReadAsArray()
        PenIMGs[:,:,id] = cmip_data
        PenMask[cmip_data == -9999] = -9999
    pendata = PenIMGs.mean(axis=2)
    pendata[PenMask == -9999] = -9999
    out_path = os.path.join(OutDir,"cmip5.{}.{}.{}.tif".format(str(yy),str(mm).zfill(2),pentandID))
    print("d",out_path)
    Write(pendata, out_path, RefProj, RefGeoTrans, RefW, RefH, im_bands=1, dtype=gdal.GDT_Float32)

for year in range(2007,2019):
    for month in range(1,13):
        for penid in range(1,7):

            GenCmip5Pentand(r"D:\Cornell\EthiopianDrought\CMIP5DailyClip\Big",r"D:\Cornell\EthiopianDrought\CMIP5Pentand",year,month,penid)

