from osgeo import gdal,osr,ogr
from GetPointFromShp import *
from SifRetrieval import *
import xlrd
from scipy import stats
path = r"D:\cropyieldforecast\boundarypolygon.shp"


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
        gmtemp =[]
        for ge in gm:
            gmtemp.append(ge.split(" "))

        info_dict["lonlat"] = gmtemp
        cropset.append(info_dict)

    return cropset

def retrieval(spec_path,crop_path):
    Cropset = getCrop(crop_path)
    Dataset,HeadTitle,White = fileload(spec_path)
    Dataset2, HeadTitle2 = fileload2(spec_path)

    extent = (HeadTitle > 759) & (HeadTitle < 764)
    ext = np.where(extent)

    Extent, Index = getMinValue(HeadTitle, White)
    lamd = HeadTitle[Extent][Index]

    mWhite = White[ext]
    mData = []
    for data in Dataset:
        data = np.array(data)
        mData.append(data[ext])
    mData = np.array(mData)
    mHT = HeadTitle[ext]
    R = [data/mWhite for data in mData]
    R = np.array(R)
    Lonyc = []
    Latyc = []
    wLat = np.where(HeadTitle2 == "S-Latitude")
    wLon = np.where(HeadTitle2 == "S-Longitude")
    for idx,r in enumerate(R):

        if mHT[np.argmin(r)] == lamd:
            print(mHT[np.argmin(r)],r)
            print("**", mHT[np.argmin(r)], lamd)

            Lonyc.append(float(Dataset2[idx][wLon]))
            Latyc.append(float(Dataset2[idx][wLat]))





    print(Latyc)
    fig = plt.figure()
    Crop_Sif = []
    assert len(Dataset) == len(Dataset2), "mispatch"
    for cps in Cropset:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for lonlat in cps["lonlat"]:
            lon,lat = float(lonlat[0]),float(lonlat[1])
            ring.AddPoint(lon,lat)

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        lon,lat = [float(ll[0]) for ll in cps["lonlat"]],[float(ll[1]) for ll in cps["lonlat"]]

        plt.plot(lon,lat)
        plt.scatter([float(d[wLon]) for d in Dataset2],[float(d[wLat]) for d in Dataset2] )
        plt.scatter(Lonyc, Latyc,marker="s")








cp_sif = retrieval(r"D:\HUBEIFIELDDATA\20180830airline1-5\15350964458442982835611128394267_20180830113249.06.xlsx",path)


plt.show()