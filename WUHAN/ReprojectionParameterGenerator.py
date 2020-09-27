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
ImgPath = r'D:\aria2\hongboshi\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753\S2A_MSIL2A_20190105T030111_N0211_R032_T50RKU_20190105T064753.SAFE\GRANULE\L2A_T50RKU_A018476_20190105T030703\IMG_DATA\R20m\T50RKU_20190105T030111_B04_20m.jp2'




def ParameterGenerator(S2TileGranule,SrcEPSG,TargetEPSG,ParOutPath):
    st_time = time.time()
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(SrcEPSG)
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(TargetEPSG)
    ct = osr.CoordinateTransformation(src_srs, target_srs)
    ct_inv = osr.CoordinateTransformation(target_srs,src_srs)

    g = gdal.Open(S2TileGranule)
    geo_t = g.GetGeoTransform()
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize
    pixel_spacing = geo_t[1]

    # 50度带的参数
    p1x, p1y = geo_t[0], geo_t[3]
    p2x, p2y = geo_t[0] + pixel_spacing * x_size, geo_t[3]
    p3x, p3y = geo_t[0] + pixel_spacing * x_size, geo_t[3] - pixel_spacing * y_size
    p4x, p4y = geo_t[0], geo_t[3] - pixel_spacing * y_size

    end_time1 = time.time()
    print("1",end_time1 - st_time)
    # 投影后顶点坐标
    rep1x, rep1y = ct.TransformPoint(p1x, p1y)[0:2]
    rep2x, rep2y = ct.TransformPoint(p2x, p2y)[0:2]
    rep3x, rep3y = ct.TransformPoint(p3x, p3y)[0:2]
    rep4x, rep4y = ct.TransformPoint(p4x, p4y)[0:2]

    minX = min(rep1x, rep4x)
    maxX = max(rep2x, rep3x)
    minY = min(rep3y, rep4y)
    maxY = max(rep1y, rep2y)

    # 49度带的对角坐标
    ulx = int(minX / pixel_spacing) * pixel_spacing
    uly = int(np.ceil((maxY) / pixel_spacing)) * pixel_spacing
    lrx = int(np.ceil((maxX) / pixel_spacing)) * pixel_spacing
    lry = int(minY / pixel_spacing) * pixel_spacing

    end_time2 = time.time()
    print("2", end_time2 - end_time1)

    # 49度带的宽高
    Width = int((lrx - ulx)/pixel_spacing)
    Heigth = int((uly - lry)/pixel_spacing)
    print("W H",Width,Heigth)


    WD1 = ulx + np.repeat(a=[np.arange(Width)], repeats=Heigth, axis=0).flatten() * pixel_spacing + 0.5 * pixel_spacing
    HD1 = uly - np.repeat(a=np.arange(Heigth), repeats=Width) * pixel_spacing - 0.5 * pixel_spacing
    Pair = np.column_stack((WD1, HD1))

    end_time3 = time.time()
    print("3", end_time3 - end_time2)

    # 49度带转50度
    Pair = ct_inv.TransformPoints(Pair)
    end_time4 = time.time()
    print("4", end_time4 - end_time3)

    Pair = np.array(Pair).reshape(Heigth,Width,3)
    PairX = Pair[:,:,0]
    PairY = Pair[:, :, 1]
    # PairX = np.array([p[0] for p in Pair]).reshape(Heigth, Width)
    # PairY = np.array([p[1] for p in Pair]).reshape(Heigth, Width)

    # print(PairX)
    mask = (PairX < p1x) | (PairX >= p2x) | (PairY > p1y) | (PairY <= p4y)

    PairX[mask] = 0
    PairY[mask] = 0
    Valid = PairX != 0


    COL = np.floor((PairX - p1x) / pixel_spacing).astype(np.int32)
    ROW = np.floor((p1y - PairY) / pixel_spacing).astype(np.int32)

    ECord = WD1.reshape(Heigth,Width)
    NCord = HD1.reshape(Heigth, Width)

    end_time5 = time.time()
    print("5", end_time5 - end_time4)

    ParaMeterX = np.zeros(shape=(y_size,x_size),dtype=np.float)
    ParaMeterY = np.zeros(shape=(y_size,x_size),dtype=np.float)

    ValidCOL = COL[Valid]
    ValidROW = ROW[Valid]

    ValidRC = ValidCOL + ValidROW*x_size
    np.put(a=ParaMeterX, ind=ValidRC,v=ECord[Valid])
    np.put(a=ParaMeterY, ind=ValidRC, v=NCord[Valid])


    print(ParaMeterY.min())
    end_time6 = time.time()
    print("6", end_time6 - end_time5)

    RepXCord = np.zeros((Heigth,Width))
    RepYCod = np.zeros((Heigth,Width))

    Xpath = os.path.join(ParOutPath,'TestRepX{}.tif'.format(pixel_spacing))
    Ypath = os.path.join(ParOutPath,'TestRepY{}.tif'.format(pixel_spacing))
    write_Img(ParaMeterX, Xpath, src_srs, geo_t,x_size,y_size,im_bands=1, dtype=gdal.GDT_UInt32)
    write_Img(ParaMeterY, Ypath, src_srs, geo_t, x_size, y_size, im_bands=1, dtype=gdal.GDT_UInt32)

    end_time7 = time.time()
    print("7", end_time7 - end_time6)
    end_time = time.time()

    print("I am done ! take time for {} s".format(end_time - st_time))

ParameterGenerator(ImgPath,32650,32649,r'C:\Users\zmhwh\Desktop\Temp')

