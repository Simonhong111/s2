from osgeo import gdal,osr,ogr
import numpy as np
import os
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


def agg(rfpath):

    crop_path = rfpath
    reference_path = r"D:\Cornell\EthiopianDrought\CMIP5Daily\cmip5_20060101.tif"
    crop = gdal.Open(crop_path)
    crop_geo = crop.GetGeoTransform()
    cropX,cropY = crop.RasterXSize,crop.RasterYSize
    crop_proj = crop.GetProjection()
    # crop_raster = crop.ReadAsArray()
    print("crop",crop_geo)

    reference = gdal.Open(reference_path)
    reference_geo = reference.GetGeoTransform()
    reference_proj = reference.GetProjection()
    reference_raster = reference.ReadAsArray()
    Width,Height = reference.RasterXSize,reference.RasterYSize

    Aggarr = np.zeros_like(reference_raster,dtype=np.float)

    print("re",reference_geo)
    print(Height)
    print(Width)
    for row in range(Height):
        for col in range(Width):

            p1x,p1y = reference_geo[0]+col*reference_geo[1],reference_geo[3]+row*reference_geo[5]
            p2x,p2y = p1x + reference_geo[1],p1y
            p3x, p3y = p2x,p2y+reference_geo[5]
            p4x, p4y = p1x,p3y

            if p1x < crop_geo[0] or p2x > crop_geo[0] + cropX*crop_geo[1] or p1y > crop_geo[3] or p3y < crop_geo[3]+cropY*crop_geo[5]:
                continue


            C1,R1 = (p1x - crop_geo[0]) / crop_geo[1],(p1y - crop_geo[3]) / crop_geo[5]
            C2, R2 = (p2x - crop_geo[0]) / crop_geo[1], (p2y - crop_geo[3]) / crop_geo[5]
            C3, R3 = (p3x - crop_geo[0]) / crop_geo[1], (p3y - crop_geo[3]) / crop_geo[5]
            C4, R4 = (p4x - crop_geo[0]) / crop_geo[1], (p4y - crop_geo[3]) / crop_geo[5]


            W = C3-C1
            H = R3 -R1

            W = int(np.round(W))
            H = int(np.round(H))

            mC1 = int(np.round(C1))
            mR1 = int(np.round(R1))

            data = crop.ReadAsArray(mC1,mR1,W,H)
            # print(C1,C3,R1,R3)
            mask = np.where(data != -9999)
            if mask[0].shape[0] < 300:
                Aggarr[row][col] = -9999
            else:
                Aggarr[row][col] = data[data != -9999].mean()

            # Aggarr[row][col] = data[data != -9999].mean()
            # print(data[data != -9999])




    path = os.path.join(r"D:\Cornell\EthiopianDrought\CHIRPSDailyAgg",'agg_'+os.path.basename(rfpath))
    print(path)
    write_Img(Aggarr, path, reference_proj, reference_geo,Width, Height,im_bands=1, dtype=gdal.GDT_Float32)
    #
    output_raster = os.path.join(r"D:\Cornell\EthiopianDrought\CHIRPSDailyAgg",os.path.basename(rfpath))
    input_shape = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"
    clipbyshp(path,output_raster,input_shape, dstNodata=-9999)

import glob

chirps = glob.glob(r"D:\Cornell\EthiopianDrought\CHIRPSDaily\chirps-v2.0.*.tif")
for chip in chirps:
    agg(chip)