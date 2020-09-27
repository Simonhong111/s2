from osgeo import gdal,osr,ogr
import numpy as np
import os
import time

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
ImgPath = r'D:\Composite\T49RGP_201901_TCI_10m.jp2'


def ReprojecByRegcover(S2TileGranule,SrcEPSG,TargetEPSG,ParaPath,ReprojectedImgPath):
    start_time = time.time()
    st_time = time.time()
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(SrcEPSG)
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(TargetEPSG)
    ct = osr.CoordinateTransformation(src_srs, target_srs)

    g = gdal.Open(S2TileGranule)
    geo_t = g.GetGeoTransform()
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize
    pixel_spacing = geo_t[1]
    p1x, p1y = geo_t[0], geo_t[3]
    p2x, p2y = geo_t[0] + pixel_spacing * x_size, geo_t[3]
    p3x, p3y = geo_t[0] + pixel_spacing * x_size, geo_t[3] - pixel_spacing * y_size
    p4x, p4y = geo_t[0], geo_t[3] - pixel_spacing * y_size

    rep1x, rep1y = ct.TransformPoint(p1x, p1y)[0:2]
    rep2x, rep2y = ct.TransformPoint(p2x, p2y)[0:2]
    rep3x, rep3y = ct.TransformPoint(p3x, p3y)[0:2]
    rep4x, rep4y = ct.TransformPoint(p4x, p4y)[0:2]

    minX = min(rep1x, rep4x)
    maxX = max(rep2x, rep3x)
    minY = min(rep3y, rep4y)
    maxY = max(rep1y, rep2y)

    ulx = int(minX / pixel_spacing) * pixel_spacing
    uly = int(np.ceil((maxY) / pixel_spacing)) * pixel_spacing
    lrx = int(np.ceil((maxX) / pixel_spacing)) * pixel_spacing
    lry = int(minY / pixel_spacing) * pixel_spacing

    Width = int((lrx - ulx) / pixel_spacing)
    Heigth = int((uly - lry) / pixel_spacing)
    Fill_value = 20000
    geotrans = [ulx,pixel_spacing,0,uly,0,-pixel_spacing]
    print(geotrans,Width,Heigth)
    t1 = time.time()
    print('time1',t1-st_time)
    RepImg = np.zeros(shape=(Heigth,Width), dtype=np.int16)
    Para = gdal.Open(ParaPath)
    t2 = time.time()
    print('time2', t2 - t1)
    PairX = Para.GetRasterBand(1).ReadAsArray()
    PairY = Para.GetRasterBand(2).ReadAsArray()
    PairXF = PairX[PairX !=Fill_value].astype(np.int64)
    PairYF = PairY[PairY !=Fill_value].astype(np.int64)

    del PairY
    Ind = PairYF * Width + PairXF

    t3 = time.time()
    print('time3', t3 - t2)
    raster = g.ReadAsArray()
    t4 = time.time()
    print('time4', t4 - t3)

    np.put(a=RepImg,ind=Ind,v=raster[:,:,0][PairX !=Fill_value])



    del PairX
    del Ind
    t5 = time.time()
    print('time5', t5 - t4)
    PairXGapFiller =  Para.GetRasterBand(3).ReadAsArray()
    PairYGapFiller =  Para.GetRasterBand(4).ReadAsArray()

    t6 = time.time()
    print('time6', t6 - t5)
    EpValidCOL = PairXGapFiller[PairXGapFiller !=Fill_value].astype(np.int64)
    EpValidROW = PairYGapFiller[PairYGapFiller != Fill_value].astype(np.int64)

    t7 = time.time()
    print('time7', t7 - t6)

    del PairYGapFiller
    EpValidRC = EpValidCOL + EpValidROW * Width
    print("E",EpValidRC.max())
    del EpValidCOL
    del EpValidROW

    np.put(RepImg, EpValidRC, raster[PairXGapFiller != Fill_value])

    del raster
    del EpValidRC
    del PairXGapFiller
    t8 = time.time()
    print('time8', t8 - t7)
    Path = os.path.join(ReprojectedImgPath,os.path.basename(S2TileGranule)[:-4]+'b2.0.tif')
    write_Img(RepImg, Path, target_srs, geotrans, Width, Heigth, im_bands=1, dtype=gdal.GDT_Float32)
    t9 = time.time()
    print('time8', t9 - t8)

    end_time = time.time()
    print("I am done and take time {} s".format(end_time - start_time))

PairPath = r'C:\Users\zmhwh\Desktop\Temp\T49RGP_TO50N_R10_E119510_N3404700_W11550_H11550.tif'

ReprojectedImgPath = r'D:\Composite'
ReprojecByRegcover(ImgPath,32649,32650,PairPath,ReprojectedImgPath)
# raster = gdal.Open(r'C:\Users\zmhwh\Desktop\Temp\T50RKU_20190105T030111_B04_60mRep.tif').ReadAsArray()
# print(raster[11,0])