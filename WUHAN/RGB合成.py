from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob
import time
from multiprocessing import Pool
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



def RGBComposite(TileName,Tm,TileDir,OutDir):
    st = time.time()
    #L2A_201912_T49RGQ_20200426T222639_R10m_NIR
    # R = glob.glob(os.path.join(TileDir,'L2A_{}_{}_*_R10m_B04.tif').format(Tm,TileName))[0]
    # G = glob.glob(os.path.join(TileDir, 'L2A_{}_{}_*_R10m_B03.tif').format(Tm,TileName))[0]
    # B = glob.glob(os.path.join(TileDir, 'L2A_{}_{}_*_R10m_B02.tif').format(Tm,TileName))[0]
    R = os.path.join(TileDir,'T50RKU_20190418T030551_B04_10m.jp2')
    G = os.path.join(TileDir, 'T50RKU_20190418T030551_B03_10m.jp2')
    B = os.path.join(TileDir, 'T50RKU_20190418T030551_B02_10m.jp2')
    assert os.path.exists(R),"no file name {} exists".format(R)
    assert os.path.exists(G), "no file name {} exists".format(G)
    assert os.path.exists(B), "no file name {} exists".format(B)


    Rr = gdal.Open(R)
    Gr = gdal.Open(G)
    Br = gdal.Open(B)
    proj = Rr.GetProjection()
    geotrans = Rr.GetGeoTransform()
    basename = os.path.basename(R).split("_")
    # Name = basename[0]+'_'+basename[1]+'_'+basename[2]+'_'+basename[3]+'_'+basename[4]+'_'+'RGB'+'.tif'
    Name = 'testf8.jpg'
    Path = os.path.join(OutDir,Name)
    print(Path)

    Size = 10980


    Data = np.zeros(shape=(Size,Size,3))
    Data[:,:,0] = Rr.ReadAsArray()*0.0255
    del Rr
    Data[:, :, 1] = Gr.ReadAsArray()*0.0255
    del Gr
    Data[:, :, 2] = Br.ReadAsArray()*0.0255
    del Br
    write_Img(Data, Path, proj, geotrans, 10980, 10980, im_bands=3, dtype=gdal.GDT_Byte)

    end = time.time()
    print(end - st)



TileName = 'T50RKU'
Tm = '201901'

TileDir = r'D:\Composite\hh'
Res = 10
OutDir = r'C:\Users\zmhwh\Desktop\Temp'
# RGBComposite(TileName,Tm,TileDir,OutDir)
path = r'C:\Users\zmhwh\Desktop\Temp\t.jp2'
path2 = r'C:\Users\zmhwh\Desktop\Temp\t2.jp2'

raster = gdal.Open(r'D:\Composite\hh\T50RKU_20190418T030551_TCI_10m.jp2')

help(raster)
print("*",raster.GetGCPSpatialRef())
print(raster.GetGCPCount())
print(raster.GetGCPProjection())
print(raster.GetGCPs())
# print(raster.GetMetadataItem())