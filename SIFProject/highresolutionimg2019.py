from osgeo import gdal,osr,ogr
from GetPointFromShp import *
from SifRetrieval2019 import *
import xlrd
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
path = r"D:\Satellive\downshp\downland.shp"

def GetFieldGeoInfo(res,spec_path):
    """
    :param croppath: 这个是实验天的路径
    :param res: 设置需要生成的影像的分辨率
    :return: 返回值Info 设计到，矢量的范围（矩形），凸包，范围投影到对应UTM坐标带的geotransform 信息，投影信息，
    """
    Dataset, HeadTitle, White = fileload(spec_path)
    Dataset2, HeadTitle2 = fileload2(spec_path)

    # 获取经纬度对应的列
    wLat = np.where(HeadTitle2 == "S-Latitude")
    wLon = np.where(HeadTitle2 == "S-Longitude")
    wAlt = np.where(HeadTitle2 == "S-Altitude(m)")
    UAVLon = [float(data[wLon[0]]) for data in Dataset2]
    UAVLat = [float(data[wLat[0]]) for data in Dataset2]

    multipoint = ogr.Geometry(ogr.wkbMultiPoint)


    for idx, lon in enumerate(UAVLon):
        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(lon, UAVLat[idx])
        multipoint.AddGeometry(point1)  #
        point1 = None

    mpConvexHull = multipoint.ConvexHull()  # 探针采样点的凸包
    gm = str(mpConvexHull)[10:-2].split(",")
    # print("gm",gm)
    Lon2 = []
    Lat2 = []
    for ge in gm:
        Lat2.append(float(ge.split(" ")[1]))
        Lon2.append(float(ge.split(" ")[0]))
    Lon2 = np.array(Lon2)
    Lat2 = np.array(Lat2)

    ulx = Lon2.min()
    uly = Lat2.max()
    lrx = Lon2.max()
    lry = Lat2.min()
    extent =[ulx,lrx,lry,uly]
    cenLon = (extent[0] + extent[1]) / 2.0
    cenLat = (extent[2] + extent[3]) / 2.0
    # print("center Lon and Lat is ",cenLon, cenLat)

    zone = np.round((183 + cenLon) / 6, 0)
    # print("zone",zone)

    EPSG = 32700 - np.round((45 + cenLat) / 90, 0) * 100 + np.round((183 + cenLon) / 6, 0)
    EPSG = int(EPSG)
    # print("EPSG",EPSG)

    prosrs = osr.SpatialReference()
    prosrs.ImportFromEPSG(EPSG)
    geosrs = prosrs.CloneGeogCS()

    ulx, uly = extent[0], extent[3]
    lrx, lry = extent[1], extent[2]

    ct = osr.CoordinateTransformation(geosrs, prosrs)
    ulx1, uly1 = ct.TransformPoint(uly, ulx)[:2]  # 坐标系转化之后的值 注，有的版本要把Easting 放在第一个参数
    lrx1, lry1 = ct.TransformPoint(lry, lrx)[:2]

    ulx2 = np.floor(ulx1).astype(np.int)
    uly2 = np.ceil(uly1).astype(np.int)

    lrx2 = np.ceil(lrx1).astype(np.int)
    lry2 = np.floor(lry1).astype(np.int)

    # print("Top E N ",ulx1,uly1)
    # print("Bottom E N ",lrx1,lry1)
    #
    # print("Adjusted Top E N ", ulx2, uly2)
    # print("Adjusted Bottom E N ", lrx2, lry2)
    #
    info = {}

    info["geotrsform"] = [ulx2, res, 0, uly2, 0, (-1.0) * res]
    info["projection"] = prosrs
    info["geocs"] = geosrs
    info["extent"] = [ulx2, lrx2, lry2, uly2]  # Easting 最小 - > 最大 North 最小 ->
    info["convexhull"] = mpConvexHull
    # print(info)
    return info

def isContained(FieldGeoInfo,res):
    """
    :param FieldGeoInfo: 获取试验地坐标信息
    :param res:
    :return:
    """
    ulx, lrx, lry, uly = FieldGeoInfo["extent"]
    H = int(np.ceil((uly - lry)/res))  # 生成影像的长宽
    W = int(np.ceil((lrx - ulx)/res))
    HresImg = np.zeros((H,W)).astype(np.float)
    geotr = FieldGeoInfo["geotrsform"]
    convexhull = FieldGeoInfo["convexhull"]
    proj = FieldGeoInfo["projection"]
    geocs = FieldGeoInfo["geocs"]
    ct = osr.CoordinateTransformation(proj, geocs)

    # 将落入实验天中的数据都设置为99，相当于生成了试验田的掩膜
    for h in range(H):
        for w in range(W):
            Easting = geotr[0] + geotr[1] * w
            Northing = geotr[3] + geotr[5] * h
            # print(Easting,Northing)

            lat,lon =ct.TransformPoint(float(Easting) ,float(Northing))[:2]
            # print("lon lat is ",lon,lat)
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(lon,lat)
            if convexhull.Contains(point):
                HresImg[h,w] = 99
            point = None
    return HresImg

def checkValidaPoints(prosrs, geosrs,convexhull,lon,lat,sif,altitude):
    """
    :param prosrs: 投影
    :param geosrs:大地坐标
    :param convexhull:试验田凸包
    :param lon:探针采样点经度位置
    :param lat:探针采样点经度位置
    :param sif:反演的sif
    :param altitude:探针高度
    :return:
    """
    Easting = []
    Northing = []
    SIF  = []
    ALT = []
    ct = osr.CoordinateTransformation(geosrs,prosrs)
    # 获取再凸包里面的数据
    for idx,mlon in enumerate(lon):
        ulx,uly = ct.TransformPoint(float(lat[idx]),float(mlon))[:2]
        # print("ulx,uly",ulx,uly)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(mlon, lat[idx])
        if convexhull.Contains(point):
            # print("ulx{},uly{} is within the extent".format(ulx,uly))
            Easting.append(ulx)
            Northing.append(uly)
            SIF.append(sif[idx])
            ALT.append(altitude[idx])
        point = None
    return np.array(Easting),np.array(Northing),np.array(SIF),np.array(ALT)



def retrieval(spec_path):
    """
    :param spec_path: 探针采集的数据
    :return:
    """
    Dataset, HeadTitle, White = fileload(spec_path)
    Dataset2, HeadTitle2 = fileload2(spec_path)

    # 获取经纬度对应的列
    wLat = np.where(HeadTitle2 == "S-Latitude")
    wLon = np.where(HeadTitle2 == "S-Longitude")
    wAlt = np.where(HeadTitle2 == "S-Altitude(m)")
    UAVLon = [float(data[wLon[0]]) for data in Dataset2]
    UAVLat = [float(data[wLat[0]]) for data in Dataset2]
    UAVAlt = [float(data[wAlt[0]]) for data in Dataset2]
    UAVLat = np.array(UAVLat).astype(np.float)
    UAVLon = np.array(UAVLon).astype(np.float)
    UAVAlt = np.array(UAVAlt).astype(np.float)
    SIF = []

    for idx, mlat in enumerate(UAVLat):
        SIF.append(Rtr3FLD_AVE(Dataset[idx], HeadTitle, White))  # FLD反演的荧光
        # print(RtrsFLD(Dataset[idx],HeadTitle,White,LeftStart=757))
    return UAVLon,UAVLat, UAVAlt, np.array(SIF).astype(np.float)

def NeighborPixel(FieldGeoInfo,East,North,ALt,SIF,res):
    """
    :param FieldGeoInfo:
    :param East:
    :param North:
    :param ALt:
    :param SIF:
    :param res:
    :return: 计算每一个采样点数据，所对应的位置
    """
    ulx, lrx, lry, uly = FieldGeoInfo["extent"]
    H = int(np.ceil((uly - lry) / res))
    W = int(np.ceil((lrx - ulx) / res))

    East2  = []
    North2 = []
    Alt2 = []
    SIF2 = []
    for idx, e in enumerate(East):
        dis = np.tan(12.5/180*np.pi)*ALt[idx] * 0.5
        LE = int(np.floor(max(0,e-dis))/res)
        RE = int(np.ceil(min(e + dis, W))/res)
        UN = int(np.floor(max(0, North[idx] - dis))/res)
        BN = int(np.ceil(min(North[idx] + dis, H))/res)
        # print(LE,RE,UN,BN)
        for w in range(LE,RE):
            for h in  range(UN,BN):
                if ((w+0.5)*res - e) * ((w+0.5)*res - e) + ((h + 0.5)*res - North[idx]) * \
                        ((h + 0.5) * res - North[idx]) < dis*dis:
                    # print("dddd")
                    East2.append(w)
                    North2.append(h)
                    Alt2.append(ALt[idx])
                    SIF2.append(SIF[idx])
    return np.array(East2),np.array(North2),np.array(Alt2),np.array(SIF2)
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


def GenImgMatrix(spec_path,res,DIS=1):

    FieldGeoInfo = GetFieldGeoInfo(res,spec_path)
    prosrs = FieldGeoInfo["projection"]
    geosrs = FieldGeoInfo["geocs"]
    geotr = FieldGeoInfo["geotrsform"]
    convexhull = FieldGeoInfo["convexhull"]
    HreImg = isContained(FieldGeoInfo,res) # 试验田掩膜
    lon,lat, altitude,sif = retrieval(spec_path) # 从探针数据获取的信息

    Easting,Northing,SIF,ALT =\
        checkValidaPoints(prosrs, geosrs, convexhull, lon, lat, sif, altitude)
    Cols = (Easting - geotr[0])/res
    Rows = (geotr[3] - Northing)/res
    East2, North2, Alt2, SIF2 = NeighborPixel(FieldGeoInfo,Easting - geotr[0],geotr[3] - Northing,ALT,SIF,res)
    # 获取了实验田的范围，获取了有效点，接下来就是插值了

    ImgSum = np.zeros((HreImg.shape[0],HreImg.shape[1])) # 把SIF写入图像，但是由于肯能有重叠，需要求均值
    ImgNum = np.zeros((HreImg.shape[0],HreImg.shape[1]))
    multipoint = ogr.Geometry(ogr.wkbMultiPoint)

    # SS = 0
    # NN = 0

    for idx, e in enumerate(East2):

        ImgSum[North2[idx],e] += SIF2[idx]
        ImgNum[North2[idx],e] += 1
        # if e ==37 and North2[idx] ==38:
        #     SS += SIF2[idx]
        #     NN += 1
        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(float(e), float(North2[idx]))
        multipoint.AddGeometry(point1)  #
        point1 = None
    # print("SS/NN",SS/NN)
    mpConvexHull = multipoint.ConvexHull()   # 探针采样点的凸包
    Mask = np.where(ImgNum ==0)
    ImgNum[Mask] = -1
    ImgSum[Mask] = 99
    # print(ImgSum)

    ResImg = ImgSum/ImgNum

    airArea = HreImg.copy()
    H,W = HreImg.shape
    E = []
    N = []
    for h in range(H):
        for w in range(W):
            tempoint = ogr.Geometry(ogr.wkbPoint)
            tempoint.AddPoint(w,h)

            if mpConvexHull.Contains(tempoint):
                ROWs = []
                COLs = []
                VALs = []
                E.append(w)
                N.append(h)
                if ImgNum[h,w] == -1:
                    for row in range(h-DIS,h+DIS+1):
                        for col in range(w-DIS,w+DIS+1):
                            if ImgNum[row,col] != -1:
                                ROWs.append(row)
                                COLs.append(col)
                                VALs.append(ResImg[row,col])
                weigthsum = 0
                val = 0

                if len(COLs) > 1:
                    ImgNum[h,w] = 1
                    for ix, mcol in enumerate(COLs):

                        val += ResImg[ROWs[ix],mcol]/(np.square(ROWs[ix] - h)+ np.square(mcol - w))
                        weigthsum += 1/(np.square(ROWs[ix] - h)+ np.square(mcol - w))

                    ResImg[h,w] = val/weigthsum
    output = spec_path[:-5] + ".tif"
    write_Img(ResImg,output,prosrs,geotr,W,H)
    # 需要一个测试的
    # fig1  = plt.figure(1)
    # plt.imshow(ImgNum)
    # fig1 = plt.figure(2)
    # fig2 = plt.imshow(ResImg)
    # # plt.scatter(Cols,Rows,c="r")
    #
    # plt.scatter(East2,North2,c="b")
    # # plt.scatter(E,N,marker='o',c='',edgecolors='k')
    # print(len(East2))
    # plt.show()







spec_dir = r"D:\Satellive\20190517-20190827\tozhang"
sub_spec = os.listdir(spec_dir)
sub_spec = [s for s in sub_spec if s.startswith("2019")]
# print(sub_spec)
p_spec = [os.path.join(spec_dir,s) for s in sub_spec]
# print(p_spec)
import  glob

for ssub in p_spec:
    # print(ssub)
    spec_fns = glob.glob(os.path.join(ssub,"A*"))
    print(spec_fns)
    for fns in spec_fns:

        spec_xlsx = glob.glob(os.path.join(fns, "*.xlsx"))
        print(spec_xlsx)
        for s in spec_xlsx:
            print("***",s)
            GenImgMatrix(s,1)



