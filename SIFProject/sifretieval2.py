from osgeo import gdal,osr,ogr
from GetPointFromShp import *
from SifRetrieval2019 import *
import xlrd
from scipy import stats
path = r"D:\Satellive\cropyieldforecast\boundarypolygon.shp"
spec_path = \
    r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180828airline1-5\15350964458442982835611128394267_20180828110010.64.xlsx"

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
def visulize(spec_path,path):
    Dataset, HeadTitle, White = fileload(spec_path)
    Dataset2, HeadTitle2 = fileload2(spec_path)

    # 获取经纬度对应的列
    wLat = np.where(HeadTitle2 == "S-Latitude")
    wLon = np.where(HeadTitle2 == "S-Longitude")
    wAlt = np.where(HeadTitle2 == "S-Altitude(m)")
    print(wLat)
    UAVLon = [float(data[wLon[0]]) for data in Dataset2]
    UAVLat = [float(data[wLat[0]]) for data in Dataset2]
    UAVAlt = [float(data[wAlt[0]]) for data in Dataset2]

    fig = plt.figure()

    plt.scatter(UAVLon,UAVLat)

    cropset = getCrop(path)
    for cps in cropset:
        # print(cps["lonlat"])
        Lon = [float(LL[0]) for LL in cps["lonlat"]]
        Lat = [float(LL[1]) for LL in cps["lonlat"]]
        plt.plot(Lon,Lat)
    plt.show()

visulize(spec_path,path)

# path = r"D:\Satellive\uppershp\upperland.shp"
# path = r"D:\Satellive\downshp\downland.shp"
# Data = getCrop(path)
#
# Yield = []
# for data in Data:
#     temp = []
#     temp.append(float(data["yield"]))
#     temp.append(float(data["area"]))
#     Yield.append(temp)
# print(Yield)
# Yield = np.array(Yield).astype(np.float)
# YieldPer = Yield[:,0]/Yield[:,1]
# print(YieldPer.mean())

# Ａ１　　２０７５　　２２４５
# Ａ２　　２２４２　　　２５３０
# Ａ３１　　１５１５　　　１７６５
# Ａ３２　　９４３　　　　９５５
# Ａ４　　　２６２０　　　３２０５
# Ａ５　　　　４２５０　　４６２０
# Ａ６　　　３４６０　　　３８３０

# Area = np.array([2075,2242,1515,943,2620,4250,3460])
# Yd =   np.array([2245,2530,1765,955,3205,4620,3830])
# YP = Yd/Area
# print(YP.mean())












