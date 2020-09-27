import os,glob
import numpy as np
from osgeo import gdal,osr,ogr
from geo2mapxy import *
from Sentinel2Reader import SenL2AReader

def SaveSen2Agg(masked_grid_path,sen2dir,save_npz_path):

    """
    :param grid_path: 读取有效的网格数据
    :param sen2dir: 哨兵影像数据文件.safe 文件
    :param save_npz_path: 存储数据对
    :return:
    """

    mGridInfo = np.load(masked_grid_path)
    print("mGridInfo",mGridInfo.files)
    mGridId = mGridInfo["sifgridId"]
    mGridTile = mGridInfo["gridTile"]
    mGridEtRing = mGridInfo["EtRing"]
    print("LENGTH",len(mGridId))


    GridID = []
    GridEtRing = []
    Date = []
    BandArr = []
    GridS2Path =[]
    gridNum = 0

    for idx,mRing in enumerate(mGridEtRing):

        if mGridTile[idx].__len__() >= 1:


            senlist = readSen2ByMonth(sen2dir,mGridTile[idx])

            if senlist is not None and senlist != []:
                print("**",senlist)
                ulx = float(mRing[0])  # 计算矢量的范围，左上角右下角
                uly = float(mRing[1])
                lrx = float(mRing[4])
                lry = float(mRing[5])

                for sen2 in senlist:
                    # print(ulx,uly)

                    mSenL2A = SenL2AReader(sen2)
                    cloudPro = mSenL2A.getQiData("MSK_CLDPRB")
                    ulxgeo,ulygeo = lonlat2geo(cloudPro,uly,ulx)
                    lrxgeo,lrygeo = lonlat2geo(cloudPro,lry,lrx)

                    ulcol,ulrow = geo2imagexy(cloudPro,ulxgeo,ulygeo)
                    lrcol, lrrow = geo2imagexy(cloudPro, lrxgeo, lrygeo)

                    width = int(np.ceil(lrcol - ulcol))
                    heigth = int(np.ceil(lrrow - ulrow))

                    # print(width, heigth)

                    cloud = cloudPro.ReadAsArray(ulcol,ulrow,width,heigth).astype(np.int)

                    cloudmask = np.where(cloud < 30)


                    if len(cloudmask[0]) > 0 and ((cloudmask[0].__len__()*1.0)/cloud.size>0.97):


                        GridID.append(mGridId[idx])
                        Date.append(os.path.basename(sen2).split("_")[2][0:8])
                        GridEtRing.append(mGridEtRing[idx])
                        GridS2Path.append(sen2)

                        bandmean = []
                        cldmsk = cloud < 30

                        bandset =["B02","B03","B04","B05","B06","B07","B8A","B11","B12"]

                        for band in bandset:
                            raster = mSenL2A.getImgData(band).ReadAsArray(ulcol,ulrow,width,heigth).astype(np.int)

                            mask = (raster < 10000) & (raster >0) & cldmsk
                            mask = np.where(mask == True)

                            # print("mask",mask)
                            if mask[0].__len__() > (raster.size*0.75) :

                                mean = np.mean(raster[mask])*0.0001
                                bandmean.append(mean)
                                # print(band, mean)

                                del raster
                                del mask
                        BandArr.append(bandmean)


                    del mSenL2A
                    del cloudPro

    kwargs = {"gridId": GridID, "gridEtRing": GridEtRing, "Date": Date,"mean":BandArr,"s2path":GridS2Path}
    np.savez(save_npz_path,**kwargs)



def readSen2ByMonth(sen2dir,tilename):
    senlist = []
    for name in tilename:
        temp = glob.glob(os.path.join(sen2dir,"T"+name,"*L2A_20180*.SAFE"))
        senlist.extend(temp)
    return senlist
