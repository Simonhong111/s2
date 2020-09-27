from osgeo import gdal,osr,ogr
from GetPointFromShp import *
from SifRetrieval import *
import xlrd
from scipy import stats
import pandas as pd
path = r"D:\Satellive\cropyieldforecast\boundarypolygon.shp"


def genImgInfo(path,res):

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(path, 0)
    layer = dataSource.GetLayer()

    extent = layer.GetExtent()
    print("img extent",extent)
    cenLon = (extent[0] + extent[1])/2.0
    cenLat = (extent[2] + extent[3])/2.0
    print("cenLon,cenLat",cenLon,cenLat)

    zone = np.round((183+cenLon)/6,0)
    print("zone",zone)
    EPSG=32700-np.round((45+cenLat)/90,0)*100+np.round((183+cenLon)/6,0)
    EPSG = int(EPSG)
    print("EPSG",EPSG)

    prosrs = osr.SpatialReference()
    prosrs.ImportFromEPSG(EPSG)
    geosrs = prosrs.CloneGeogCS()

    ulx,uly = extent[0],extent[3]
    lrx,lry = extent[1],extent[2]

    ct = osr.CoordinateTransformation(geosrs,prosrs)
    ulx1,uly1 = ct.TransformPoint(uly,ulx)[:2]
    lrx1,lry1 = ct.TransformPoint(lry,lrx)[:2]

    print(ulx1, uly1)
    print(lrx1, lry1)
    ulx2, uly2 = np.floor(ulx1) - 10, np.ceil(uly1) + 10
    lrx2, lry2 = np.ceil(lrx1) + 10 , np.floor(lry1) - 10
    print(ulx2, uly2)
    print(lrx2, lry2)

    nCols = int(np.ceil((lrx2-ulx2)/res))
    nRows = int(np.ceil((uly2-lry2)/res))
    info ={}

    info["geotrsform"] = [ulx2,res,0,uly2,0,(-1.0)*res]
    info["projection"] = prosrs
    info["geocs"] = geosrs
    info["ncols"] = nCols
    info["nrows"] = nRows

    return info

def getCrop(crop_shp_path):

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(path, 0)
    layer = dataSource.GetLayer()
    cropid = 0
    cropset =[]
    for feature in layer:
        cropid += 1
        info_dict = {}
        info_dict["cropId"] =  cropid
        info_dict["yield"] = feature.GetField("yield")
        info_dict['area'] = feature.GetField("area")
        temp =[]
        geom = feature.GetGeometryRef()
        # print(geom)
        gm =  str(geom)[10:-2].split(",")
        # print("gm",gm)
        gmtemp =[]
        for ge in gm:
            gmtemp.append(ge.split(" "))

        info_dict["lonlat"] = gmtemp
        cropset.append(info_dict)

    return cropset


def genImg(shp_path,res,measurement,output):

    mInfo = genImgInfo(shp_path,res)

    sif = np.zeros([mInfo["nrows"],mInfo["ncols"]])
    nMatrix = np.zeros([mInfo["nrows"],mInfo["ncols"]])

    ulx,uly = mInfo["geotrsform"][0],mInfo["geotrsform"][3]
    print("ulx",ulx,uly)

    lrx = ulx + mInfo["geotrsform"][1] * mInfo["ncols"]
    lry = uly + mInfo["geotrsform"][5] * mInfo["nrows"]

    dataset, HeadTitle, White, Dark = fileload(measurement)

    dataset = np.array(dataset)
    HT = np.array(HeadTitle)

    bandIndex = np.where(HT == '769.739')
    print("bandIndx",bandIndex)

    latIndx = np.where(HT == "S-Latitude")
    lonIndx = np.where(HT == "S-Longitude")

    lat = [float(dat[latIndx]) for dat in dataset]
    lon = [float(dat[lonIndx]) for dat in dataset]
    band = [float(dat[bandIndex]) for dat in dataset]

    lat = np.array(lat)
    lon = np.array(lon)
    band = np.array(band)

    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(str(mInfo["projection"]))
    geosrs = prosrs.CloneGeogCS()
    ct = osr.CoordinateTransformation(geosrs, prosrs)

    for idx,mlon in enumerate(lon):
        print("**",lat[idx],mlon)
        mX, mY = ct.TransformPoint(lat[idx], mlon)[:2]
        print("mx",type(mX),mX,mY)
        print("ii",ulx,uly)
        print("ii", lrx, lry)
        print(lat[idx], mlon)
        print("ture",(mX > ulx and mX < lrx))
        if (mX > ulx and mX < lrx) and (mY > lry and mY < uly):
            print("mmm",mX,mY)

            icol = int((mX - ulx)/res)
            irow = int((uly - mY)/res)

            sif[irow,icol] +=  band[idx]

            nMatrix[irow,icol] += 1

    mask = np.where(nMatrix < 1)

    sif[mask] = -1

    nMatrix[mask] = 1

    msif = sif/nMatrix


    write_Img(msif,output,mInfo["projection"],mInfo["geotrsform"],im_width=mInfo["ncols"],im_heigth=mInfo["nrows"])




def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):
    from osgeo import osr
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype)

    dataset.SetGeoTransform(geotrans)
    print("proj",proj)
    dataset.SetProjection(str(proj))
    if (im_bands == 1) :
        print("**********")
        dataset.GetRasterBand(1).WriteArray(data)
    del dataset




def retrieval(spec_path,crop_path):
    Cropset = getCrop(crop_path)
    Dataset,HeadTitle,White = fileload(spec_path)

    Dataset2, HeadTitle2 = fileload2(spec_path) # 包含了所有的列数据

    # 获取经纬度对应的列
    wLat = np.where(HeadTitle2=="S-Latitude")
    wLon = np.where(HeadTitle2=="S-Longitude")

    # fig = plt.figure()
    Crop_Sif = []
    assert len(Dataset) == len(Dataset2), "mispatch occurs"
    for cps in Cropset:#读取每一个田块的信息包括产量，位置范围

        ring = ogr.Geometry(ogr.wkbLinearRing)
        for lonlat in cps["lonlat"]:
            lon,lat = float(lonlat[0]),float(lonlat[1])
            ring.AddPoint(lon,lat)#将田块的多边形点存入环中

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)#将环存入多边形中

        obser_num = 0
        Sum_FLD, Sum_FLD_AVE,Sum_3FLD,Sum_3FLD_AVE = 0.0,0.0,0.0,0.0
        mSum_FLD, mSum_FLD_AVE, mSum_3FLD, mSum_3FLD_AVE = -99, -99, -99, -99

        Sum_NDVI,Sum_EVI,Sum_GRVI,Sum_GNDVI,Sum_MSAVI,Sum_OSAVI,Sum_SAVI,Sum_WDRVI = 0,0,0,0,0,0,0,0
        mSum_NDVI,mSum_EVI,mSum_GRVI,mSum_GNDVI,mSum_MSAVI,mSum_OSAVI,mSum_SAVI,mSum_WDRVI = \
            -99,-99,-99,-99,-99,-99,-99,-99   # 这些是为了求均值的呢

        for idx in range(len(Dataset2)):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(float(Dataset2[idx][wLon]),float(Dataset2[idx][wLat]))

            if poly.Contains(point):
                print("Contains",cps["cropId"],cps["yield"],cps["area"],cps["lonlat"],Dataset[idx][0])

                Sum_FLD += RtrsFLD(Dataset[idx],HeadTitle,White,LeftStart=757) # FLD反演的荧光

                Sum_FLD_AVE += RtrsFLD_AVE(Dataset[idx],HeadTitle,White) # FLD平滑后反演的荧光

                Sum_3FLD += Rtr3FLD(Dataset[idx],HeadTitle,White) # 3FLD反演的荧光

                Sum_3FLD_AVE += Rtr3FLD_AVE(Dataset[idx], HeadTitle, White) # 3FLD平滑后反演的荧光
                NDVI2, EVI2, GRVI2, GNDVI2, MSAVI2, OSAVI2, SAVI2, WDRVI2 = VIS(Dataset[idx], HeadTitle)

                # Sum_NDVI += NDVI(Dataset[idx], HeadTitle)  # NDVI
                # Sum_EVI      += EVI(Dataset[idx], HeadTitle)
                # Sum_GRVI     += GRVI(Dataset[idx], HeadTitle)
                # Sum_GNDVI    += GNDVI(Dataset[idx], HeadTitle)
                # Sum_MSAVI    += MSAVI(Dataset[idx], HeadTitle)
                # Sum_OSAVI    += OSAVI(Dataset[idx], HeadTitle)
                # Sum_SAVI     += SAVI(Dataset[idx], HeadTitle)
                # Sum_WDRVI    += WDRVI(Dataset[idx], HeadTitle)

                Sum_NDVI += NDVI2
                Sum_EVI      += EVI2
                Sum_GRVI     += GRVI2
                Sum_GNDVI    += GNDVI2
                Sum_MSAVI    += MSAVI2
                Sum_OSAVI    += OSAVI2
                Sum_SAVI     += SAVI2
                Sum_WDRVI    += WDRVI2

                obser_num  += 1  #观测数量+1

        # obser_num = 0
        # Sum_FLD1, Sum_FLD_AVE1,Sum_3FLD1,Sum_3FLD_AVE1,Sum_NDVI1 = [],[],[],[],[]
        # mSum_FLD, mSum_FLD_AVE, mSum_3FLD, mSum_3FLD_AVE = -99, -99, -99, -99
        # mSum_NDVI = -99   # 这些是为了求均值的呢
        #
        # for idx in range(len(Dataset2)):
        #     point = ogr.Geometry(ogr.wkbPoint)
        #     point.AddPoint(float(Dataset2[idx][wLon]),float(Dataset2[idx][wLat]))
        #
        #     if poly.Contains(point):
        #         print("Contains",cps["cropId"],Dataset[idx][0])
        #
        #         Sum_FLD1.append(RtrsFLD(Dataset[idx],HeadTitle,White,LeftStart=757)) # FLD反演的荧光
        #
        #         Sum_FLD_AVE1.append(RtrsFLD_AVE(Dataset[idx],HeadTitle,White))  # FLD平滑后反演的荧光
        #
        #         Sum_3FLD1.append(Rtr3FLD(Dataset[idx],HeadTitle,White) )# 3FLD反演的荧光
        #
        #         Sum_3FLD_AVE1.append(Rtr3FLD_AVE(Dataset[idx], HeadTitle, White))  # 3FLD平滑后反演的荧光
        #
        #         Sum_NDVI1.append(NDVI(Dataset[idx], HeadTitle))   # NDVI
        #
        #         obser_num  += 1  #观测数量+1


        if obser_num < 5:
            mSum_FLD = -99
            mSum_FLD_AVE = -99
            mSum_3FLD = -99
            mSum_3FLD_AVE = -99
            mSum_NDVI = -99
            mSum_EVI  = -99
            mSum_GRVI = -99
            mSum_GNDVI = -99
            mSum_MSAVI = -99
            mSum_OSAVI = -99
            mSum_SAVI = -99
            mSum_WDRVI = -99

        else:
            # Sum_FLD =     np.array(Sum_FLD1).sum() - np.array(Sum_FLD1).min() - np.array(Sum_FLD1).max()
            # Sum_FLD_AVE = np.array(Sum_FLD_AVE1).sum() - np.array(Sum_FLD_AVE1).min() - np.array(Sum_FLD_AVE1).max()
            # Sum_3FLD =    np.array(Sum_3FLD1).sum() - np.array(Sum_3FLD1).min() - np.array(Sum_3FLD1).max()
            # Sum_3FLD_AVE = np.array(Sum_3FLD_AVE1).sum() - np.array(Sum_3FLD_AVE1).min()  -np.array(Sum_3FLD_AVE1).max()
            # Sum_NDVI = np.array(Sum_NDVI1).sum() - np.array(Sum_NDVI1).min() - np.array(Sum_NDVI1).max()
            # obser_num = obser_num - 2
            mSum_FLD = Sum_FLD/obser_num   # mSum_FLD means the average of the sum
            mSum_FLD_AVE = Sum_FLD_AVE/obser_num
            mSum_3FLD = Sum_3FLD/obser_num
            mSum_3FLD_AVE = Sum_3FLD_AVE/obser_num
            mSum_NDVI = Sum_NDVI/obser_num
            mSum_EVI = Sum_EVI/obser_num
            mSum_GRVI = Sum_GRVI /obser_num
            mSum_GNDVI = Sum_GNDVI/obser_num
            mSum_MSAVI = Sum_MSAVI/obser_num
            mSum_OSAVI = Sum_OSAVI/obser_num
            mSum_SAVI = Sum_SAVI/obser_num
            mSum_WDRVI = Sum_WDRVI/obser_num

        info_dict ={}
        info_dict["cropId"] = cps["cropId"]
        info_dict["yield"] = cps["yield"]
        info_dict['area'] = cps['area']
        info_dict["FLD"] = mSum_FLD
        info_dict["FLD_AVE"] = mSum_FLD_AVE
        info_dict["3FLD"] = mSum_3FLD
        info_dict["3FLD_AVE"] = mSum_3FLD_AVE
        info_dict["NDVI"] = mSum_NDVI
        info_dict["EVI"] = mSum_EVI
        info_dict["GRVI"] = mSum_GRVI
        info_dict["GNDVI"] = mSum_GNDVI
        info_dict["MSAVI"] = mSum_MSAVI
        info_dict["OSAVI"] = mSum_OSAVI
        info_dict["SAVI"] = mSum_SAVI
        info_dict["WDRVI"] = mSum_WDRVI

        PP = str(poly.Centroid())[7:-2].split(" ")

        info_dict["Lon"] = float(PP[0])
        info_dict["Lat"] =float(PP[1])

        Crop_Sif.append(info_dict)

        # mLon,mLat = [float(ll[0]) for ll in cps["lonlat"]],[float(ll[1]) for ll in cps["lonlat"]]
        # plt.plot(mLon,mLat)
        # plt.text(float(PP[0]), float(PP[1]), 'sif = {}'.format(round(mSum_FLD,3)))
        # plt.scatter([float(d[wLon]) for d in Dataset2],[float(d[wLat]) for d in Dataset2] )
    print(Crop_Sif)
    return Crop_Sif








def plotD(measurements,path):

    cp_sif = retrieval(measurements, path)
    Area = np.array([float(cps['area']) for cps in cp_sif])
    FLD = np.array([float(cps["FLD"]) for cps in cp_sif])
    FLD_AVE = np.array([float(cps["FLD_AVE"]) for cps in cp_sif])
    TrFLD = np.array([float(cps["3FLD"]) for cps in cp_sif])
    TrFLD_AVE = np.array([float(cps["3FLD_AVE"]) for cps in cp_sif])
    NDVI = np.array([float(cps["NDVI"]) for cps in cp_sif])
    EVI = np.array([float(cps["EVI"]) for cps in cp_sif])
    GRVI = np.array([float(cps["GRVI"]) for cps in cp_sif])
    GNDVI = np.array([float(cps["GNDVI"]) for cps in cp_sif])
    MSAVI = np.array([float(cps["MSAVI"]) for cps in cp_sif])
    OSAVI = np.array([float(cps["OSAVI"]) for cps in cp_sif])
    SAVI = np.array([float(cps["SAVI"]) for cps in cp_sif])
    WDRVI = np.array([float(cps["WDRVI"]) for cps in cp_sif])

    YLD = np.array([float(cps["yield"]) for cps in cp_sif])

    masks = np.where(FLD != -99)
    mFLD, mFLD_AVE, m3FLD, m3FLD_AVE,mNDVI,mEVI,mGRVI ,mGNDVI ,mMSAVI ,mOSAVI ,mSAVI ,mWDRVI = FLD[masks], FLD_AVE[masks], TrFLD[masks], TrFLD_AVE[masks],NDVI[masks], \
                                             EVI[masks],GRVI[masks] ,GNDVI[masks] ,MSAVI[masks] ,OSAVI[masks] ,SAVI[masks] ,WDRVI[masks]
    Area, YLD = Area[masks], YLD[masks]

    mFLD = mFLD * Area * 0.001
    mFLD_AVE = mFLD_AVE * Area * 0.001
    m3FLD = m3FLD * Area * 0.001
    m3FLD_AVE = m3FLD_AVE * Area * 0.001
    mNDVI = mNDVI*Area*0.001
    mEVI  = mEVI*Area*0.001
    mGRVI = mGRVI*Area*0.001
    mGNDVI = mGNDVI*Area*0.001
    mMSAVI = mMSAVI*Area*0.001
    mOSAVI = mOSAVI*Area*0.001
    mSAVI = mSAVI*Area*0.001
    mWDRVI = mWDRVI*Area*0.001


    print("NDVI",mNDVI)

    slope1, intercept1, r_value_fld, p_value_fld, std_err = stats.linregress(YLD, mFLD)
    print("FLD cor", r_value_fld ** 2,slope1)

    slope2, intercept2, r_value_fld2, p_value_fld2, std_err = stats.linregress(YLD, mFLD_AVE)
    print("FLD_AVE cor", r_value_fld2 ** 2,slope2)

    slope3, intercept3, r_value_3fld, p_value_3fld, std_err = stats.linregress(YLD, m3FLD)
    print("3FLD cor", r_value_3fld ** 2,slope3)

    slope4, intercept4, r_value_3fld2, p_value_3fld2, std_err = stats.linregress(YLD, m3FLD_AVE)
    print("3FLD_AVE cor", r_value_3fld2 ** 2,slope4)

    slope5, intercept5, r_value_ndvi, p_value_ndvi, std_err = stats.linregress(YLD, mNDVI)
    print("NDVI cor", r_value_ndvi ** 2,slope5)

    slope6, intercept6, r_value_evi, p_value_evi, std_err = stats.linregress(YLD, mEVI)
    print("EVI cor", r_value_evi ** 2, slope6)

    slope7, intercept7, r_value_grvi, p_value_grvi, std_err = stats.linregress(YLD, mGRVI)
    print("GRVI cor", r_value_grvi ** 2, slope7)

    slope8, intercept8, r_value_gndvi, p_value_gndvi, std_err = stats.linregress(YLD, mGNDVI)
    print("GNDVI cor", r_value_gndvi ** 2, slope8)

    slope9, intercept9, r_value_msavi, p_value_msavi, std_err = stats.linregress(YLD, mMSAVI)
    print("MSAVI cor", r_value_msavi ** 2, slope9)

    slope10, intercept10, r_value_osavi, p_value_osavi, std_err = stats.linregress(YLD, mOSAVI)
    print("OSAVI cor", r_value_osavi ** 2, slope10)

    slope11, intercept11, r_value_savi, p_value_savi, std_err = stats.linregress(YLD, mSAVI)
    print("SAVI cor", r_value_savi ** 2, slope11)

    slope12, intercept12, r_value_wdrvi, p_value_wdrvi, std_err = stats.linregress(YLD, mWDRVI)
    print("WDRVI cor", r_value_wdrvi ** 2, slope12)



    xmin,xmax = np.min(YLD) ,np.max(YLD)
    print("dat",os.path.basename(measurements).split("_")[1][:-5])
    month = os.path.basename(measurements).split("_")[1][:-5][5:6]
    day = os.path.basename(measurements).split("_")[1][:-5][6:8]


    if day[0] =="0":
        day = day[1]
    print(month,day)
    # dataframe = pd.DataFrame(
    #     {"date":str(os.path.basename(measurements).split("_")[1][:-5]),"CropArea": Area, 'Yield': YLD, "FLD": mFLD, "FLD_AVE": mFLD_AVE, "3FLD": m3FLD, "3FLD_AVE": m3FLD_AVE,"NDVI":mNDVI,\
    #      "r2_1":r_value_fld, "p1":p_value_fld,"s1":slope1,"inc1":intercept1,\
    #      "r2_2":r_value_fld2, "p2":p_value_fld2,"s2":slope2,"inc2":intercept2,\
    #      "r2_3":r_value_3fld, "p3":p_value_3fld,"s3":slope3,"inc3":intercept3,\
    #      "r2_4":r_value_3fld2, "p4":p_value_3fld2,"s4":slope4,"inc4":intercept4,\
    #      "r2_5":r_value_ndvi, "p5":p_value_ndvi,"s5":slope5,"inc5":intercept5,})

    # 将DataFrame存储为csv,index表示是否显示行名，default=True
    # dataframe.to_csv(os.path.join(r"D:\HUBEIFIELDDATA\Results", "maifei.csv"), index=False, sep=',',mode="a")
    fig = plt.figure(figsize=(20,8))
    fig.add_subplot(3,4,1)
    plt.title("FLD R2={}, p={}".format(round(r_value_fld ** 2, 3), round(p_value_fld, 3)))

    plt.scatter(YLD, mFLD)
    plt.plot([xmin,xmax],[slope1*xmin + intercept1,slope1*xmax + intercept1],"r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope1 + intercept1, "{}*x+{}".format(round(slope1,3), round(intercept1,3)))

    fig.add_subplot(3,4,2)
    plt.title("FLD_AVE r2={},p={}".format(round(r_value_fld2 ** 2, 3), round(p_value_fld2, 3)))
    plt.scatter(YLD, mFLD_AVE)
    plt.plot([xmin,xmax],[slope2*xmin + intercept2,slope2*xmax + intercept2],"r")
    plt.text(xmin/2+xmax/2,(xmin/2+xmax/2)*slope2 + intercept2,"{}*x+{}".format(round(slope2,3), round(intercept2,3)))

    fig.add_subplot(3,4,3)
    plt.title("3FLD r2={},p={}".format(round(r_value_3fld ** 2, 3), round(p_value_3fld, 3)))
    plt.scatter(YLD, m3FLD)
    plt.plot([xmin,xmax],[slope3*xmin + intercept3,slope3*xmax + intercept3],"r")
    plt.text(xmin/2+xmax/2,(xmin/2+xmax/2)*slope3 + intercept3,"{}*x+{}".format(round(slope3,3), round(intercept3,3)))

    fig.add_subplot(3,4,4)
    plt.title("3FLD_AVE r2={},p={}".format(round(r_value_3fld2 ** 2, 3), round(p_value_3fld2, 3)))
    plt.scatter(YLD, m3FLD_AVE)
    plt.plot([xmin,xmax],[slope4*xmin + intercept4,slope4*xmax + intercept4],"r")
    plt.text(xmin/2+xmax/2,(xmin/2+xmax/2)*slope4 + intercept4,"{}*x+{}".format(round(slope4,3), round(intercept4,3)))

    fig.add_subplot(3,4,5)
    plt.title("NDVI r2={},p={}".format(round(r_value_ndvi ** 2, 3), round(p_value_ndvi, 3)))
    plt.scatter(YLD, mNDVI)
    plt.plot([xmin,xmax],[slope5*xmin + intercept5,slope5*xmax + intercept5],"r")
    plt.text(xmin/2+xmax/2,(xmin/2+xmax/2)*slope5 + intercept5,"{}*x+{}".format(round(slope5,3), round(intercept5,3)))

    fig.add_subplot(3,4,6)
    plt.title("EVI r2={},p={}".format(round(r_value_evi ** 2, 3), round(p_value_evi, 3)))
    plt.scatter(YLD, mEVI)
    plt.plot([xmin, xmax], [slope6 * xmin + intercept6, slope6 * xmax + intercept6], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope6 + intercept6,
             "{}*x+{}".format(round(slope6, 3), round(intercept6, 3)))

    fig.add_subplot(3,4,7)
    plt.title("GRVI r2={},p={}".format(round(r_value_grvi ** 2, 3), round(p_value_grvi, 3)))
    plt.scatter(YLD, mGRVI)
    plt.plot([xmin, xmax], [slope7 * xmin + intercept7, slope7 * xmax + intercept7], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope7 + intercept7,
             "{}*x+{}".format(round(slope7, 3), round(intercept7, 3)))

    fig.add_subplot(3,4,8)
    plt.title("GNDVI r2={},p={}".format(round(r_value_gndvi ** 2, 3), round(p_value_gndvi, 3)))
    plt.scatter(YLD, mGNDVI)
    plt.plot([xmin, xmax], [slope8 * xmin + intercept8, slope8 * xmax + intercept8], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope8 + intercept8,
             "{}*x+{}".format(round(slope8, 3), round(intercept8, 3)))




    fig.add_subplot(3,4,9)
    plt.title("MSAVI r2={},p={}".format(round(r_value_msavi ** 2, 3), round(p_value_msavi, 3)))
    plt.scatter(YLD, mMSAVI)
    plt.plot([xmin, xmax], [slope9 * xmin + intercept9, slope9 * xmax + intercept9], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope9 + intercept9,
             "{}*x+{}".format(round(slope9, 3), round(intercept9, 3)))


    fig.add_subplot(3,4,10)
    plt.title("OSAI r2={},p={}".format(round(r_value_osavi ** 2, 3), round(p_value_osavi, 3)))
    plt.scatter(YLD, mOSAVI)
    plt.plot([xmin, xmax], [slope10 * xmin + intercept10, slope10 * xmax + intercept10], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope10 + intercept10,
             "{}*x+{}".format(round(slope10, 3), round(intercept10, 3)))


    fig.add_subplot(3,4,11)
    plt.title("SAVI r2={},p={}".format(round(r_value_savi ** 2, 3), round(p_value_savi, 3)))
    plt.scatter(YLD, mSAVI)
    plt.plot([xmin, xmax], [slope11 * xmin + intercept11, slope11 * xmax + intercept11], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope11 + intercept11,
             "{}*x+{}".format(round(slope11, 3), round(intercept11, 3)))

    fig.add_subplot(3,4,12)
    plt.title("WDRVI r2={},p={}".format(round(r_value_wdrvi ** 2, 3), round(p_value_wdrvi, 3)))
    plt.scatter(YLD, mWDRVI)
    plt.plot([xmin, xmax], [slope12 * xmin + intercept12, slope12 * xmax + intercept12], "r")
    plt.text(xmin / 2 + xmax / 2, (xmin / 2 + xmax / 2) * slope12 + intercept12,
             "{}*x+{}".format(round(slope12, 3), round(intercept12, 3)))




    fig.savefig(os.path.join(os.path.dirname(measurements), os.path.basename(measurements).split("_")[1][:-5] + ".png"))

    # plt.plot(np.arange(len(mYld)),mYld,c="b",label="CropYield")
    # plt.plot(np.arange(len(mYld)),mSIF*5,c="k",label="SIF_FLD")
    # plt.plot(np.arange(len(mYld)),mSIF_3FLD*5,c="g",label="SIF_3FLD")
    # plt.plot(np.arange(len(mYld)),mNDVI*500,c="r",label="NDVI")

    plt.legend()
    plt.show()

import glob

subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180825airline6-10","*.xlsx"))
for d in subdir:

    plotD(d,path)

# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180828airline1-5","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
#
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180830airline1-5","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180830airline6-10","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180903airline1-5","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180905airline6-10","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180914airline6-10","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180915airline6-10","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
# subdir = glob.glob(os.path.join(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180918airline1-5","*.xlsx"))
# for d in subdir:
#
#     plotD(d,path)
