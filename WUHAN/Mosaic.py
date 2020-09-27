from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob
from matplotlib import pyplot as plt


def getextent(boundaryshp,baseepsg,pixel_spacing):
    DriverName = "ESRI Shapefile"
    driver = ogr.GetDriverByName(DriverName)
    if driver is None:
        print("%s driver not available.\n" % boundaryshp)
        return
    else:
        print("%s driver IS available.\n" % boundaryshp)

    dataSource = driver.Open(boundaryshp, 0)
    layer = dataSource.GetLayer()
    extent = layer.GetExtent() # [minlon,maxlon,minlat,maxlat]

    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    basesrs = osr.SpatialReference()
    basesrs.ImportFromEPSG(baseepsg)
    ct = osr.CoordinateTransformation(wgs84,basesrs)

    LatLon =[[extent[3],extent[0]],[extent[3],extent[1]],[extent[2],extent[1]],[extent[2],extent[0]]]
    LatLon =  np.array(LatLon)
    PairsRep = ct.TransformPoints(LatLon)
    PairsRep = np.array(PairsRep).reshape(len(LatLon),3)
    Ecord,Ncord = PairsRep[:,0],PairsRep[:,1]
    minX,maxX,minY,maxY =  Ecord.min(),Ecord.max(),Ncord.min(),Ncord.max()
    # 49度带的对角坐标
    ulx = int(minX / pixel_spacing) * pixel_spacing
    uly = int(np.ceil((maxY) / pixel_spacing)) * pixel_spacing
    lrx = int(np.ceil((maxX) / pixel_spacing)) * pixel_spacing
    lry = int(minY / pixel_spacing) * pixel_spacing
    Width = int((lrx - ulx) / pixel_spacing)
    Heigth = int((uly - lry) / pixel_spacing)

    return [ulx,uly,lrx,lry,Width,Heigth,pixel_spacing]
def Insect(T1Name,T2Name,T1,T2):
    # Create ring
    ring1 = ogr.Geometry(ogr.wkbLinearRing)
    ring1.AddPoint(T1['P1X'],T1['P1Y'])
    ring1.AddPoint(T1['P2X'], T1['P2Y'])
    ring1.AddPoint(T1['P3X'], T1['P3Y'])
    ring1.AddPoint(T1['P4X'], T1['P4Y'])
    ring1.AddPoint(T1['P1X'], T1['P1Y'])
    # Create polygon
    poly1 = ogr.Geometry(ogr.wkbPolygon)
    poly1.AddGeometry(ring1)

    ring2 = ogr.Geometry(ogr.wkbLinearRing)
    ring2.AddPoint(T2['P1X'],T2['P1Y'])
    ring2.AddPoint(T2['P2X'], T2['P2Y'])
    ring2.AddPoint(T2['P3X'], T2['P3Y'])
    ring2.AddPoint(T2['P4X'], T2['P4Y'])
    ring2.AddPoint(T2['P1X'], T2['P1Y'])
    # Create polygon
    poly2 = ogr.Geometry(ogr.wkbPolygon)
    poly2.AddGeometry(ring2)

    if poly1.Intersect(poly2):
        print("{} and {} Intesect".format(T1Name,T2Name))



def segment(TileDir,Date,Paras,Res=10,Size = 10980):

    T50RKV = glob.glob(os.path.join(TileDir, 'T50RKV_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T50RLV = glob.glob(os.path.join(TileDir, 'T50RLV_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T50RKU = glob.glob(os.path.join(TileDir, 'T50RKU_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T49RGP = glob.glob(os.path.join(TileDir, 'T49RGP_{}*_TCI_{}m.tif'.format(Date, Res)))[0]

    T50RKV = gdal.Open(T50RKV)
    T50RLV = gdal.Open(T50RLV)
    T50RKU = gdal.Open(T50RKU)


    T50RKVGeot = T50RKV.GetGeoTransform()
    T50RLVGeot = T50RLV.GetGeoTransform()
    T50RKUGeot = T50RKU.GetGeoTransform()

    T50RKVP1X, T50RKVP1Y = T50RKVGeot[0], T50RKVGeot[3]
    T50RKVP2X, T50RKVP2Y = T50RKVP1X + Res * Size, T50RKVP1Y
    T50RKVP3X, T50RKVP3Y = T50RKVP2X, T50RKVP2Y - Res * Size
    T50RKVP4X, T50RKVP4Y = T50RKVP1X, T50RKVP3Y

    T50RLVP1X, T50RLVP1Y = T50RLVGeot[0], T50RLVGeot[3]
    T50RLVP2X, T50RLVP2Y = T50RLVP1X + Res * Size, T50RLVP1Y
    T50RLVP3X, T50RLVP3Y = T50RLVP2X, T50RLVP2Y - Res * Size
    T50RLVP4X, T50RLVP4Y = T50RLVP1X, T50RLVP3Y

    T50RKUP1X, T50RKUP1Y = T50RKUGeot[0], T50RKUGeot[3]
    T50RKUP2X, T50RKUP2Y = T50RKUP1X + Res * Size, T50RKUP1Y
    T50RKUP3X, T50RKUP3Y = T50RKUP2X, T50RKUP2Y - Res * Size
    T50RKUP4X, T50RKUP4Y = T50RKUP1X, T50RKUP3Y

    W = 11550
    H = 11550
    T49RGPP1X, T49RGPP1Y = 119510, 3404700
    T49RGPP2X, T49RGPP2Y = T49RGPP1X + Res * W, T49RGPP1Y
    T49RGPP3X, T49RGPP3Y = T49RGPP2X, T49RGPP2Y - Res * H
    T49RGPP4X, T49RGPP4Y = T49RGPP1X, T49RGPP3Y



    TileDict = {}
    TileDict["T50RKV"] = {"P1X":T50RKVP1X,  "P1Y":T50RKVP1Y,"P2X":T50RKVP2X,  "P2Y":T50RKVP2Y,
                          "P3X":T50RKVP3X,  "P3Y":T50RKVP3Y,"P4X":T50RKVP4X,  "P4Y":T50RKVP4Y,}
    TileDict["T50RLV"] = {"P1X": T50RLVP1X, "P1Y": T50RLVP1Y, "P2X": T50RLVP2X, "P2Y": T50RLVP2Y,
                          "P3X": T50RLVP3X, "P3Y": T50RLVP3Y, "P4X": T50RLVP4X, "P4Y": T50RLVP4Y, }
    TileDict["T50RKU"] = {"P1X": T50RKUP1X, "P1Y": T50RKUP1Y, "P2X": T50RKUP2X, "P2Y": T50RKUP2Y,
                          "P3X": T50RKUP3X, "P3Y": T50RKUP3Y, "P4X": T50RKUP4X, "P4Y": T50RKUP4Y, }

    TileDict["T49RGP"] = {"P1X":T49RGPP1X , "P1Y": T49RGPP1Y, "P2X": T49RGPP2X, "P2Y": T49RGPP2Y,
                          "P3X": T49RGPP3X, "P3Y": T49RGPP3Y, "P4X": T49RGPP4X, "P4Y": T49RGPP4Y, }

    ulx, uly, lrx, lry, Width, Heigth, pixel_spacing = Paras
    print(lrx)
    TileDict['T50RLV']['P1X'] = int((TileDict['T50RKV']['P2X'] + TileDict['T50RLV']['P1X'])/2) -2
    TileDict['T50RLV']['P1Y'] = uly + 1
    TileDict['T50RLV']['P2Y'] = uly + 1
    TileDict['T50RLV']['P4X'] = TileDict['T50RLV']['P1X']
    TileDict['T50RLV']['P2X'] = lrx +2
    TileDict['T50RLV']['P3X'] = lrx + 2

    TileDict['T50RKV']['P1Y'] = uly + 2
    TileDict['T50RKV']['P2Y'] = uly + 2
    TileDict['T50RKV']['P2X'] = TileDict['T50RLV']['P1X'] +2
    TileDict['T50RKV']['P3X'] = TileDict['T50RLV']['P1X'] +2
    TileDict['T50RKV']['P3Y'] = int((TileDict['T50RKV']['P3Y'] + TileDict['T50RKU']['P1Y'])/2) -4
    TileDict['T50RKV']['P4Y'] = TileDict['T50RKV']['P3Y']-4


    TileDict['T50RKU']['P1X'] =  int((TileDict['T50RKU']['P1X']+TileDict['T49RGP']['P2X'])/2) -2
    TileDict['T50RKU']['P4X'] =  TileDict['T50RKU']['P1X']
    TileDict['T50RKU']['P1Y'] = TileDict['T50RKV']['P4Y'] + 2
    TileDict['T50RKU']['P2Y'] = TileDict['T50RKV']['P4Y'] + 2
    TileDict['T50RKU']['P3Y'] = lry -2
    TileDict['T50RKU']['P4Y'] = lry -2

    TileDict['T49RGP']['P1X'] = ulx -2
    TileDict['T49RGP']['P4X'] = ulx - 2
    TileDict['T49RGP']['P2X'] = TileDict['T50RKU']['P1X'] +2
    TileDict['T49RGP']['P3X'] = TileDict['T50RKU']['P1X'] +2
    TileDict['T49RGP']['P3Y'] = lry -2
    TileDict['T49RGP']['P4Y'] = lry -2


    whX,whY =[ulx,lrx,lrx,ulx,ulx],[uly,uly,lry,lry,uly]
    T50RLVX,T50RLVY = [TileDict['T50RLV']['P1X'],TileDict['T50RLV']['P2X'],TileDict['T50RLV']['P3X'],TileDict['T50RLV']['P4X'],
                       TileDict['T50RLV']['P1X']], [TileDict['T50RLV']['P1Y'],TileDict['T50RLV']['P2Y'],TileDict['T50RLV']['P3Y'],TileDict['T50RLV']['P4Y'],
                       TileDict['T50RLV']['P1Y']]

    T50RKVX, T50RKVY = [TileDict['T50RKV']['P1X'], TileDict['T50RKV']['P2X'], TileDict['T50RKV']['P3X'],
                        TileDict['T50RKV']['P4X'],
                        TileDict['T50RKV']['P1X']], [TileDict['T50RKV']['P1Y'], TileDict['T50RKV']['P2Y'],
                                                     TileDict['T50RKV']['P3Y'], TileDict['T50RKV']['P4Y'],
                                                     TileDict['T50RKV']['P1Y']]

    T50RKUX, T50RKUY = [TileDict['T50RKU']['P1X'], TileDict['T50RKU']['P2X'], TileDict['T50RKU']['P3X'],
                        TileDict['T50RKU']['P4X'],
                        TileDict['T50RKU']['P1X']], [TileDict['T50RKU']['P1Y'], TileDict['T50RKU']['P2Y'],
                                                     TileDict['T50RKU']['P3Y'], TileDict['T50RKU']['P4Y'],
                                                     TileDict['T50RKU']['P1Y']]

    T49RGPX, T49RGPY = [TileDict['T49RGP']['P1X'], TileDict['T49RGP']['P2X'], TileDict['T49RGP']['P3X'],
                        TileDict['T49RGP']['P4X'],
                        TileDict['T49RGP']['P1X']], [TileDict['T49RGP']['P1Y'], TileDict['T49RGP']['P2Y'],
                                                     TileDict['T49RGP']['P3Y'], TileDict['T49RGP']['P4Y'],
                                                     TileDict['T49RGP']['P1Y']]
    # plt.plot(whX,whY)
    # plt.plot(T50RLVX,T50RLVY)
    # plt.plot(T50RKVX, T50RKVY)
    # plt.plot(T50RKUX, T50RKUY)
    # plt.plot(T49RGPX, T49RGPY)
    # plt.show()

    return  TileDict["T50RKV"],TileDict["T50RLV"],TileDict["T50RKU"], TileDict['T49RGP']

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


def mosaic(TileDir,OutputDir,Wuhanshp,Date,Res=10,Size = 10980,OutputName =None):
    # L2A_201901_T49RGQ_20200426T210807_R10m_B02

    Paras = getextent(Wuhanshp, 32650, 10)
    T50RKV,T50RLV,T50RKU, T49RGP = segment(TileDir, Date,Paras, Res=10, Size=10980)

    T50RKVP = glob.glob(os.path.join(TileDir, 'T50RKV_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T50RLVP = glob.glob(os.path.join(TileDir, 'T50RLV_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T50RKUP = glob.glob(os.path.join(TileDir, 'T50RKU_{}*_TCI_{}m.jp2'.format(Date, Res)))[0]
    T49RGPP = glob.glob(os.path.join(TileDir, 'T49RGP_{}*_TCI_{}m.tif'.format(Date, Res)))[0]

    T50RKVR = gdal.Open(T50RKVP)
    T50RLVR = gdal.Open(T50RLVP)
    T50RKUR = gdal.Open(T50RKUP)
    T49RGPR = gdal.Open(T49RGPP)

    [ulx, uly, lrx, lry, Width, Heigth, pixel_spacing] = Paras

    ulx = ulx - 10*pixel_spacing
    uly = uly + 10*pixel_spacing
    Width = Width + 50
    Heigth = Heigth + 50
    geo_t = [ulx,pixel_spacing,0,uly,0,-pixel_spacing]
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(32650)
    Img = np.zeros(shape=(Heigth,Width,3),dtype=np.float)
    print(Heigth,Width)
    # 写入49RGP
    mW,mH = int((T49RGP["P2X"]-T49RGP["P1X"])/pixel_spacing),int((T49RGP["P1Y"]-T49RGP["P4Y"])/pixel_spacing)
    geot49 =  T49RGPR.GetGeoTransform()
    ulx49 = geot49[0]
    uly49 = geot49[3]
    col = int((T49RGP["P1X"] - ulx49)/pixel_spacing)
    row = int((uly49 - T49RGP["P1Y"])/pixel_spacing)
    T49RGPR = T49RGPR.ReadAsArray(col,row,mW,mH)
    mcol = int((T49RGP["P1X"] - ulx)/pixel_spacing)
    mrow = int((uly - T49RGP["P1Y"]) / pixel_spacing)
    Img[mrow:mrow+mH,mcol:mcol+mW,0] = T49RGPR[0,:,:]
    Img[mrow:mrow + mH, mcol:mcol + mW, 1] = T49RGPR[1, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 2] = T49RGPR[2, :, :]

    del T49RGPR

    # 写入50RKU
    mW, mH = int((T50RKU["P2X"] - T50RKU["P1X"]) / pixel_spacing), int((T50RKU["P1Y"] - T50RKU["P4Y"]) / pixel_spacing)
    geot49 = T50RKUR.GetGeoTransform()
    ulx49 = geot49[0]
    uly49 = geot49[3]
    col = int((T50RKU["P1X"] - ulx49) / pixel_spacing)
    row = int((uly49 - T50RKU["P1Y"]) / pixel_spacing)
    T50RKUR = T50RKUR.ReadAsArray(col, row, mW, mH)
    mcol = int((T50RKU["P1X"] - ulx) / pixel_spacing)
    mrow = int((uly - T50RKU["P1Y"]) / pixel_spacing)
    Img[mrow:mrow + mH, mcol:mcol + mW, 0] = T50RKUR[0, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 1] = T50RKUR[1, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 2] = T50RKUR[2, :, :]

    del T50RKUR

    # 写入49RGP
    mW, mH = int((T50RKV["P2X"] - T50RKV["P1X"]) / pixel_spacing), int((T50RKV["P1Y"] - T50RKV["P4Y"]) / pixel_spacing)
    geot49 = T50RKVR.GetGeoTransform()
    ulx49 = geot49[0]
    uly49 = geot49[3]
    col = int((T50RKV["P1X"] - ulx49) / pixel_spacing)
    row = int((uly49 - T50RKV["P1Y"]) / pixel_spacing)
    T50RKVR = T50RKVR.ReadAsArray(col, row, mW, mH)
    mcol = int((T50RKV["P1X"] - ulx) / pixel_spacing)
    mrow = int((uly - T50RKV["P1Y"]) / pixel_spacing)
    Img[mrow:mrow + mH, mcol:mcol + mW, 0] = T50RKVR[0, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 1] = T50RKVR[1, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 2] = T50RKVR[2, :, :]

    del T50RKVR

    # 写入49RGP
    mW, mH = int((T50RLV["P2X"] - T50RLV["P1X"]) / pixel_spacing), int((T50RLV["P1Y"] - T50RLV["P4Y"]) / pixel_spacing)
    geot49 = T50RLVR.GetGeoTransform()
    ulx49 = geot49[0]
    uly49 = geot49[3]
    col = int((T50RLV["P1X"] - ulx49) / pixel_spacing)
    row = int((uly49 - T50RLV["P1Y"]) / pixel_spacing)
    T50RLVR = T50RLVR.ReadAsArray(col, row, mW, mH)
    mcol = int((T50RLV["P1X"] - ulx) / pixel_spacing)
    mrow = int((uly - T50RLV["P1Y"]) / pixel_spacing)
    Img[mrow:mrow + mH, mcol:mcol + mW, 0] = T50RLVR[0, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 1] = T50RLVR[1, :, :]
    Img[mrow:mrow + mH, mcol:mcol + mW, 2] = T50RLVR[2, :, :]

    del T50RLVR

    Path = os.path.join(OutputDir,OutputName)
    write_Img(Img, Path, proj, geo_t, Width, Heigth, im_bands=3, dtype=gdal.GDT_UInt16)


TileDir = r'D:\Composite'
OutputDir = r'C:\Users\zmhwh\Desktop\Temp'
Date = '201910'
OutputName  = 'WuHanCity'+Date+'.tif'
wuhanshp = r'C:\Users\zmhwh\Downloads\CHN_adm\wuhancity.shp'
mosaic(TileDir,OutputDir,wuhanshp,Date,Res=10,Size=10980,OutputName =OutputName)


