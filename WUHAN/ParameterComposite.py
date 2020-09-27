from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob
import time
from multiprocessing import Pool
def write_Img(data, path, im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)




    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])
    del dataset


def ParameterComposite(TileName,ToTileName,TileDir,Res,OutDir):
    ParaMeterX = glob.glob(os.path.join(TileDir,'{}_TO{}_R{}_ParaX_*.tif').format(TileName,ToTileName,Res))[0]
    ParaMeterY = glob.glob(os.path.join(TileDir, '{}_TO{}_R{}_ParaY_*.tif').format(TileName, ToTileName, Res))[0]
    ParaMeterXGp = glob.glob(os.path.join(TileDir, '{}_TO{}_R{}_ParaXGap_*.tif').format(TileName, ToTileName, Res))[0]
    ParaMeterYGp = glob.glob(os.path.join(TileDir, '{}_TO{}_R{}_ParaYGap_*.tif').format(TileName, ToTileName, Res))[0]
    assert os.path.exists(ParaMeterX),"no file name {} exists".format(ParaMeterX)
    assert os.path.exists(ParaMeterY), "no file name {} exists".format(ParaMeterY)
    assert os.path.exists(ParaMeterXGp), "no file name {} exists".format(ParaMeterXGp)
    assert os.path.exists(ParaMeterYGp), "no file name {} exists".format(ParaMeterYGp)

    Xr = gdal.Open(ParaMeterX)
    Yr = gdal.Open(ParaMeterY)
    Xgpr = gdal.Open(ParaMeterXGp)
    Ygpr = gdal.Open(ParaMeterYGp)
    # T49RGP_TO50N_R10_ParaYGap_E1275660_N3432210_W11638_H11639.tif
    basename = os.path.basename(ParaMeterX).split("_")
    Name = basename[0]+'_'+basename[1]+'_'+basename[2]+'_'+basename[4]+'_'+basename[5]+'_'+basename[6]+'_'+basename[7]
    Path = os.path.join(OutDir,Name)
    print(Path)
    if Res == 10:
        Size = 10980
    elif Res == 20:
        Size = 5490
    elif Res == 60:
        Size = 1830

    Data = np.zeros(shape=(Size,Size,4))
    Data[:,:,0] = Xr.ReadAsArray()
    del Xr
    Data[:, :, 1] = Yr.ReadAsArray()
    del Yr
    Data[:, :, 2] = Xgpr.ReadAsArray()
    del Xgpr
    Data[:, :, 3] = Ygpr.ReadAsArray()
    del Ygpr
    write_Img(Data,Path,Size,Size,im_bands=4,dtype=gdal.GDT_UInt16)

TileName = 'T49RGP'
ToTileName = '50N'

TileDir = r'C:\Users\zmhwh\Desktop\Temp'
Res = 10
OutDir = r'C:\Users\zmhwh\Desktop\Temp'
ParameterComposite(TileName,ToTileName,TileDir,Res,OutDir)

