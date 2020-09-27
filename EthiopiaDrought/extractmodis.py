from osgeo import gdal,osr,ogr
from datetime import *
from dateutil import rrule
import os

# import datetime
import glob
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):
    """
    :param data: the matrix you want to save as a geotif image
    :param path: the absolute path of your image
    :param proj: the projection of your image
    :param geotrans: six parameters
    :param im_width: image width
    :param im_heigth: image height
    :param im_bands: band number of your image
    :param dtype: data type of your image
    :return:
    """

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

def dayofyear(y,m,d):
    """
    :param y: year
    :param m: month
    :param d: day
    :return: the nth day in the year
    """
    days_in_the_year = (date(y, m, d) - date(y, 1, 1)).days + 1
    return str(y)+str(days_in_the_year).zfill(3)

def mosaic(directory,outdir,VI = 1):
    """
    :param directory: where your downloaded data is stored
    :param outdir: the folder where you want to store the mosaic data
    :param VI: the index of the layers of modis product, for example ,the first layer is NDVI ,
    with VI = 0, the second layer is EVI with VI =1, the third layer is VI quality with VI = 2
    :return:
    """

    start = datetime.strptime("2000-02-01", "%Y-%m-%d").date() # time ranges of selected data that you want to process
    stop = datetime.strptime("2000-02-01", "%Y-%m-%d").date() #


    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):

        mdayofyear = dayofyear(dt.year,dt.month,dt.day) # convert your date to day of year,for example, 20030101 to 2003001

        h21v07 = glob.glob(os.path.join(directory,'*A{}.h21v07*.hdf'.format(mdayofyear)))
        h21v08 = glob.glob(os.path.join(directory,'*A{}.h21v08*.hdf'.format(mdayofyear)))
        h22v07 = glob.glob(os.path.join(directory, '*A{}.h22v07*.hdf'.format(mdayofyear)))
        h22v08 = glob.glob(os.path.join(directory,'*A{}.h22v08*.hdf'.format(mdayofyear)))

        assert len(h21v07)*len(h21v08)*len(h22v08)*len(h22v07) >0,"no such file {}".format(dt)
        ph21v07 = gdal.Open(h21v07[0]).GetSubDatasets()[VI][0]
        ph21v08 = gdal.Open(h21v08[0]).GetSubDatasets()[VI][0]
        ph22v07 = gdal.Open(h22v07[0]).GetSubDatasets()[VI][0]
        ph22v08 = gdal.Open(h22v08[0]).GetSubDatasets()[VI][0]

        print("Note***: the layer name, please check it,\n",ph21v07)
        h21v07VI = gdal.Open(ph21v07).ReadAsArray()
        h21v08VI = gdal.Open(ph21v08).ReadAsArray()
        h22v07VI = gdal.Open(ph22v07).ReadAsArray()
        h22v08VI = gdal.Open(ph22v08).ReadAsArray()

        MosaicImg = np.full(shape=(2400,2400),fill_value=-3000)  # note the fill_value is the default value of the layer ndvi and evi
        MosaicImg[0:1200,0:1200] = h21v07VI
        MosaicImg[1200:2400,0:1200] = h21v08VI
        MosaicImg[0:1200, 1200:2400] = h22v07VI
        MosaicImg[1200:2400,1200:2400] = h22v08VI

        preference = gdal.Open(ph21v07)
        proj = osr.SpatialReference()
        proj.ImportFromWkt(str(preference.GetProjection()))
        geotrans = preference.GetGeoTransform()
        path = os.path.join(outdir,str(dt.year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)+'.tif')

        write_Img(MosaicImg, path, proj, geotrans, 2400, 2400, im_bands=1, dtype=gdal.GDT_Float32)
        del MosaicImg
        del h22v08VI
        del h22v07VI
        del h21v08VI
        del h21v07VI


# mosaic(r'D:\Cornell\MOD13A3V006',r'D:\Cornell\MOD13A3V006')
# print("I am done")


def shp():
    input_shape = r'D:\Cornell\Liuyanyan\WaterShed\watshed.shp'  # this is the shapefile of watershed, see it in the folder

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(input_shape, 0)
    layer = dataSource.GetLayer()
    featureCollection = []
    # put a feature in the shapefile into dict then to a list
    for feature in layer:
        ft = {}
        # print("geo",feature.GetGeometryRef())
        # record the feature attributes and geometry
        ft['polygon'] = str(feature.GetGeometryRef())
        ft['OBJECTID'] = feature.GetField("OBJECTID")
        ft['cws_nm'] = feature.GetField("cws_nm")
        ft['mws_nm'] = feature.GetField("mws_nm")
        featureCollection.append(ft)
    return featureCollection

def point(polygon,ct):
    # get all the points in the polygon
    if 'MULTIPOLYGON' in polygon:
        polycollection = polygon.split("(((")[1].split(")))")[0].split(")),((")
        pcollection = []
        for mpy in polycollection:
            pcollection.extend(mpy.split(','))
        pairs = [[float(p.split(' ')[1]),float(p.split(' ')[0])] for p in pcollection]

    else:
        mutpoint = polygon.split("(")[2].split(")")[0].split(",")
        pairs = [[float(p.split(' ')[1]),float(p.split(' ')[0])] for p in mutpoint]
    # convert the point coordinates to a given projection
    prjpoints = ct.TransformPoints(pairs)

    East = [p[0] for p in prjpoints]
    North = [p[1] for p in prjpoints]
    return np.array(East),np.array(North)
def getRowColbyPoly():
    #
    featureCollection = shp()
    reference = gdal.Open(r"D:\Cornell\MOD13A3\20000201.tif") # this reference is the mosaic image generated by mosaic funtion mentioned above
    geo = reference.GetGeoTransform() # get the six parameters
    modissrc = reference.GetProjection() # get the projection
    targetsrc = osr.SpatialReference() # the target projection that you want to transform your data to,here it refers to modis projection
    targetsrc.ImportFromWkt(str(modissrc))
    src = osr.SpatialReference()  # the shapefile projection ,here it refers to wgs84
    src.ImportFromEPSG(4326)
    ct = osr.CoordinateTransformation(src,targetsrc) # transform coordinates in src to coordinates in targetsrc
    ct2 = osr.CoordinateTransformation(targetsrc,src)

    RC = []
    for feature in featureCollection:
        X,Y = point(feature['polygon'],ct)
        # transform the coordinates of the points of the polygon in src to the coordinates in targetsrc

        polygon = ogr.CreateGeometryFromWkt(feature['polygon'])
        # get the extent of the enclosing rectangle of points (X,Y)
        minx,maxx = X.min(),X.max()
        miny,maxy = Y.min(),Y.max()
        # convert the extent coordinates to the column and row in the mosaic coordinates
        mincol = int((minx - geo[0])/geo[1])-2
        maxcol = int((maxx - geo[0])/geo[1])+2
        minrow = int((maxy - geo[3])/geo[5])-2
        maxrow = int((miny - geo[3])/geo[5])+2

        ftRC = {}
        rlist = []
        collist = []
        lonlist = []
        latlist =[]

        # traverse every pixel in the enclosing rectangle in modis mosaic,
        # get its coordinate of its four vertices, and transform their coordinates from targetsrc(modis) to src(shapefile)
        # the use the new points to generate a polygon and decide whether the polygon of the pixel intersect the feature polygon
        # if yes, this pixel is what you want, save it column and row.
        for r in range(minrow,maxrow):
            for col in range(mincol,maxcol):
                # e - east coordinate ,n -north coordinate
                e0 = geo[0] + col*geo[1]
                n0 = geo[3] + r*geo[5]
                e1 = geo[0] + (col+1)*geo[1]
                n1 = n0
                e2 = e1
                n2 = geo[3] + (r+1)*geo[5]
                e3= e0
                n3 = n2
                lat0,lon0 = ct2.TransformPoint(e0,n0)[0:2]
                lat1,lon1 = ct2.TransformPoint(e1, n1)[0:2]
                lat2,lon2 = ct2.TransformPoint(e2, n2)[0:2]
                lat3, lon3 = ct2.TransformPoint(e3, n3)[0:2]
                wkt = 'POLYGON (({} {},{} {},{} {},{} {},{} {}))'.format(lon0,lat0,lon1,lat1,lon2,lat2,lon3,lat3,lon0,lat0)
                mpoly = ogr.CreateGeometryFromWkt(wkt)

                if polygon.Intersect(mpoly):
                    temp = []

                    rlist.append(r)
                    collist.append(col)
                    colcenter = geo[0] + col * geo[1] + geo[1]/2
                    rcenter = geo[3] + r * geo[5] + geo[5]/2
                    latC, lonC = ct2.TransformPoint(colcenter, rcenter)[0:2]
                    # covert the coordinate of center of the pixel to coordinate in src(longitude and latitude)
                    lonlist.append(lonC)
                    latlist.append(latC)


        ftRC['row'] = rlist
        ftRC['col'] = collist
        ftRC['lon']= lonlist
        ftRC['lat'] = latlist
        ftRC['polygon'] = feature['polygon']
        ftRC['OBJECTID'] = feature['OBJECTID']
        ftRC['cws_nm'] = feature['cws_nm']
        ftRC['mws_nm'] = feature['mws_nm']
        RC.append(ftRC)
    return RC

# for validation, I add the draw function for visualization
def draw(i):
    """
    :param i: i refers to which feature do you want to see
    :return:
    """
    RC = getRowColbyPoly()
    print(RC[i])

    reference = gdal.Open(r"D:\Cornell\MOD13A3\20000201.tif")
    geo = reference.GetGeoTransform()
    modissrc = reference.GetProjection()
    targetsrc = osr.SpatialReference()
    targetsrc.ImportFromWkt(str(modissrc))
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    ct = osr.CoordinateTransformation(src, targetsrc)
    ct2 = osr.CoordinateTransformation(targetsrc, src)
    east,north = point(RC[i]['polygon'],ct) # coordinate of polygon in targetsrc projection
    cx,cy = ct.TransformPoint(12.7708333321878,37.6779091378331)[0:2]  # give a point (latitide ,lontitude)
    print('cx',cx,cy)
    # convert to column and row
    cx = (cx - geo[0]) / geo[1]
    cy = (cy - geo[3]) / geo[5]
    col = (east - geo[0])/geo[1]
    row = (north - geo[3])/geo[5]
    rp = [[p,p,p+1,p+1,p] for p in RC[i]['row']] # the row of the points you get using the function of getRowColbyPoly
    cp = [[p,p+1,p+1,p,p] for p in RC[i]['col']] # the column of the points
    im = reference.ReadAsArray()
    plt.imshow(im)
    plt.scatter(cx,cy,c='y',marker='*')
    plt.plot(col,row,'r')
    for idx,r in enumerate(rp):
        plt.plot(cp[idx],r,'k')
    plt.show()

draw(12)


def extraction(directory,outpath):
    """
    :param directory: the folder where you store your  mosaic image
    :param outpath: the absolute path you want to store the csv
    :return:
    """
    start = datetime.strptime("2000-02-01", "%Y-%m-%d").date()
    stop = datetime.strptime("2000-02-01", "%Y-%m-%d").date()
    RC = getRowColbyPoly()
    PixelValues = []
    mRow= []
    mCol =[]
    mLat = []
    mLon = []
    OBJECTID = []
    cws_nm = []
    mws_nm = []
    mTime = []
    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        modis = os.path.join(directory,str(dt.year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)+'.tif')
        raster = gdal.Open(modis).ReadAsArray()
        print(os.path.basename(modis))
        for rc in RC:
            position = (np.array(rc['row']),np.array(rc['col']))
            values = raster[position]

            PixelValues.extend(values.tolist())
            mRow.extend(rc['row'])
            mCol.extend(rc['col'])
            mLat.extend(rc['lat'])
            mLon.extend(rc['lon'])
            cws_nm.extend([rc['cws_nm']]*len(rc['row']))
            mws_nm.extend([rc['mws_nm']]*len(rc['row']))
            OBJECTID.extend([rc['OBJECTID']]*len(rc['row']))
            mTime.extend([str(dt.year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)]*len(rc['row']))
    del raster
    dict = {"OBJECTID": OBJECTID, "mws_nm": mws_nm, "cws_nm":cws_nm,'row':mRow,'col':mCol,'lon':mLon,'lat':mLat,'evi':PixelValues,'time':mTime}
    eviPd = pd.DataFrame(dict)
    eviPd.to_csv(outpath)

# extraction(r'D:\Cornell\MOD13A3V006',r'D:\Cornell\MOD13A3V006\evi.csv')






























