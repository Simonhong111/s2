from osgeo import gdal,osr,ogr
import time

def GTiffTest(s2path):
    raster = gdal.Open(s2path)
    proj = raster.GetProjection()
    width,height = raster.RasterXSize,raster.RasterYSize
    bandNum = raster.RasterCount
    dtype = gdal.GDT_UInt16
    geotrans = raster.GetGeoTransform()
    path = r"D:\DriverTest\geotif_DEFLATE92.tif"
    st = time.time()
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, width, height, bandNum, dtype,options=["COMPRESS=DEFLATE"])
    dataset.SetGeoTransform(geotrans)
    dataset.SetProjection(str(proj))
    data = raster.ReadAsArray()
    if bandNum == 1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(bandNum):
            dataset.GetRasterBand(id + 1).WriteArray(data[:, :, id])
    del dataset
    et = time.time()
    print("gtiff take time for ",et-st)

def JPEG2000Test(s2path):
    raster = gdal.Open(s2path)
    proj = raster.GetProjection()
    width,height = raster.RasterXSize,raster.RasterYSize
    bandNum = raster.RasterCount
    dtype = gdal.GDT_UInt16
    geotrans = raster.GetGeoTransform()

    st = time.time()
    driver = gdal.GetDriverByName("MEM")
    dataset = driver.Create("", width, height, bandNum, dtype)
    dataset.SetGeoTransform(geotrans)
    dataset.SetProjection(str(proj))
    data = raster.ReadAsArray()
    if bandNum == 1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(bandNum):
            dataset.GetRasterBand(id + 1).WriteArray(data[:, :, id])
    JPEG_Driver = gdal.GetDriverByName("ECW")
    print("jepg_driver", JPEG_Driver)
    path = r"D:\DriverTest\JpegdriverECW .jp2"
    ds = JPEG_Driver.CreateCopy(path,dataset,0)

    del dataset
    et = time.time()
    print("gtiff take time for ",et-st)

if __name__ == '__main__':
    s2path =r'D:\Composite\hh\t49RGQ_20191111t025941_B02_10m_mosaic.tif'
    # s2path = r'D:\Composite\hh\T50RKU_20190418T030551_B04_10m.jp2'
    # s2path =r"D:\Composite\新建文件夹\S2B_MSIL2A_20191106T025909_N0213_R032_T50RKU_20191106T080154\S2B_MSIL2A_20191106T025909_N0213_R032_T50RKU_20191106T080154.SAFE\GRANULE\L2A_T50RKU_A013929_20191106T025912\IMG_DATA\R10m\T50RKU_20191106T025909_B02_10m.jp2"
    # GTiffTest(s2path)
    JPEG2000Test(s2path)



