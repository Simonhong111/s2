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
import shapefile as shp
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

def boundary(geo):
    sp = r"D:\Cornell\EthiopianDrought\ET_Region\et_region.shp"
    sf = shp.Reader(sp)
    X = []
    Y = []
    for shape in sf.shapeRecords():
        for i in range(len(shape.shape.parts)):
            i_start = shape.shape.parts[i]
            if i == len(shape.shape.parts) - 1:
                i_end = len(shape.shape.points)
            else:
                i_end = shape.shape.parts[i + 1]
            x = [(i[0]-geo[0])/geo[1] for i in shape.shape.points[i_start:i_end]]
            y = [(i[1]-geo[3])/geo[5] for i in shape.shape.points[i_start:i_end]]
            X.append(x)
            Y.append(y)
    return X,Y


def monthComposite(RainDir,RainOutDir,year,monthtype):


    if monthtype =="short":
        First = os.path.join(RainDir, "chirps-v2.0.{}.02.tif".format(str(year)))
        Second = os.path.join(RainDir, "chirps-v2.0.{}.03.tif".format(str(year)))
        Third = os.path.join(RainDir, "chirps-v2.0.{}.04.tif".format(str(year)))
        Forth = os.path.join(RainDir, "chirps-v2.0.{}.05.tif".format(str(year)))
    if monthtype == "long":
        First = os.path.join(RainDir, "chirps-v2.0.{}.06.tif".format(str(year)))
        Second = os.path.join(RainDir, "chirps-v2.0.{}.07.tif".format(str(year)))
        Third = os.path.join(RainDir, "chirps-v2.0.{}.08.tif".format(str(year)))
        Forth = os.path.join(RainDir, "chirps-v2.0.{}.09.tif".format(str(year)))

    composite = np.zeros((228, 299, 4),dtype=np.float)
    maskarr = np.zeros((228, 299))

    Img1 = gdal.Open(First).ReadAsArray()
    maskarr[Img1 == -9999] = -9999
    composite[:,:,0] = Img1
    del Img1
    Img2  = gdal.Open(Second).ReadAsArray()
    maskarr[Img2 == -9999] = -9999
    composite[:, :, 1] = Img2
    del Img2
    Img3  = gdal.Open(Third).ReadAsArray()
    maskarr[Img3 == -9999] = -9999
    composite[:, :, 2] = Img3
    del Img3
    Img4  = gdal.Open(Forth).ReadAsArray()
    maskarr[Img4 == -9999] = -9999
    composite[:, :, 3] = Img4
    del Img4

    monthComp = composite.mean(axis=2)
    monthComp[maskarr == -9999] = -9999

    chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.2010.04.tif")
    geotrans = chirps_reference.GetGeoTransform()
    proj = chirps_reference.GetProjection()
    outputpath = os.path.join(RainOutDir, monthtype + "RF" + str(year) + ".tif")
    write_Img(monthComp, outputpath, proj, geotrans, 299, 228, im_bands=1, dtype=gdal.GDT_Float32)



# for year in range(2003,2019):
#     monthComposite(r"D:\Cornell\EthiopianDrought\Chirps2",
#                      r"D:\Cornell\EthiopianDrought\AData\RainMonth",year,"long")


def countDrought(RainDir,RainOutDir,monthtype):

    ImgList = np.zeros((228, 299, 16), dtype=np.float)
    maskarr = np.zeros((228, 299))
    Drought = np.zeros((228, 299))
    for year in range(2003,2019):
        imgpath = os.path.join(RainDir,"{}RF{}.tif".format(monthtype,str(year)))
        image = gdal.Open(imgpath).ReadAsArray()
        ImgList[:,:,year-2003] = image
        maskarr[image == -9999] = -9999
        del image
    Average = ImgList.mean(axis=2)
    Std = ImgList.std(axis=2)


    for i in range(16):
        ImgList[:,:,i][maskarr != -9999] = \
            (ImgList[:,:,i][maskarr != -9999] - Average[maskarr != -9999])/Std[maskarr != -9999]
        Drought[ImgList[:,:,i] <= -1] += 1



    Drought[Drought <=1] = 1
    Drought[Drought >=4] = 4


    chirps_reference = gdal.Open(r"D:\Cornell\EthiopianDrought\Chirps2\chirps-v2.0.2010.04.tif")
    geotrans = chirps_reference.GetGeoTransform()
    proj = chirps_reference.GetProjection()
    outputpath = os.path.join(r"D:\Cornell\EthiopianDrought\AData\RainMonth", monthtype + "Drought" + str(year) + ".tif")
    write_Img(Drought, outputpath, proj, geotrans, 299, 228, im_bands=1, dtype=gdal.GDT_Float32)



    Drought[maskarr == -9999] = np.nan
    plt.title("2003-2018 Drought Frequency for Short Rains with anom < -1.0\n")
    cax = plt.imshow(Drought,cmap=plt.cm.get_cmap('coolwarm',4 ),vmin=1,vmax=5)
    # mask = np.where(Drought == 1)
    # x = mask[1]
    # y=mask[0]
    # print(mask)
    # plt.scatter(x,y,s=1,c="r")
    cbar =plt.colorbar(cax,ticks=np.arange(1,5)+0.5,label="Drought Frequency")


    cbar.ax.set_yticklabels(['<=1', '2','3','>=4'])
    plt.show()



# countDrought( r"D:\Cornell\EthiopianDrought\AData\RainMonth","","long")

def drawDrought(RainDir,RainOutDir,monthtype,myear):

    ImgList = np.zeros((228, 299, 16), dtype=np.float)
    maskarr = np.zeros((228, 299))

    for year in range(2003,2019):
        imgpath = os.path.join(RainDir,"{}RF{}.tif".format(monthtype,str(year)))
        image = gdal.Open(imgpath).ReadAsArray()
        ImgList[:,:,year-2003] = image
        maskarr[image == -9999] = -9999
        del image
    Average = ImgList.mean(axis=2)
    Std = ImgList.std(axis=2)


    for i in range(16):
        ImgList[:,:,i][maskarr != -9999] = \
            (ImgList[:,:,i][maskarr != -9999] - Average[maskarr != -9999])/Std[maskarr != -9999]


    id = myear - 2003
    Drought = ImgList[:,:,id]
    drought_id0 = Drought > -0.5
    drought_id05 = (Drought > -1) & (Drought <= -0.5)
    drought_id1 = (Drought <=-1) & (Drought >-2)
    drought_id2 = Drought <= -2

    Drought[drought_id0] = 0
    Drought[drought_id05] = -1
    Drought[drought_id1] = -2
    Drought[drought_id2] = -3
    Drought[maskarr == -9999] = np.nan



    plt.title("{} Drought Frequency for {} Rains with anom <-2,-1,-0.5\n".format(myear,monthtype))
    cax = plt.imshow(Drought,cmap=plt.cm.get_cmap('RdBu',4 ),vmin=-3,vmax=1)
    cbar =plt.colorbar(cax,ticks=np.arange(-3,1)+0.5,label="Anomaly")


    cbar.ax.set_yticklabels(['<=-2', '<=-1','<=-0.5','>-0.5'])

    geo = gdal.Open(r"D:\Cornell\EthiopianDrought\AData\RainMonth\longRF2003.tif").GetGeoTransform()
    X, Y = boundary(geo)

    for i in range(len(X)):
        plt.plot(X[i], Y[i],"r")

    # plt.tight_layout()  # 调整整体空白
    plt.show()



drawDrought( r"D:\Cornell\EthiopianDrought\AData\RainMonth","","long",2018)




