from osgeo import gdal,osr,ogr
import numpy as np
import os
import glob
import time
from multiprocessing import Pool
import shutil
import zipfile
from PIL import ImageOps,Image
def write_Img(data, path, proj, geotrans,im_width, im_heigth,im_bands=1, dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName("GTiff")

    dataset = driver.Create(path, im_width, im_heigth, im_bands, dtype,options=["COMPRESS=LZW"])

    dataset.SetGeoTransform(geotrans)
    dataset.SetProjection(str(proj))


    if im_bands ==1:
        dataset.GetRasterBand(1).WriteArray(data)
    else:
        for id in range(im_bands):
            # print("**********")
            dataset.GetRasterBand(id+1).WriteArray(data[:,:,id])

    del dataset


def decompress(zip_files, output_dir=os.getcwd(), remove=False):
    '''decompress(zip_files, output_dir = os.getcwd(), remove = False

    Decompresses .zip files downloaded from SciHub, and optionally removes original .zip file.

    Args:
        zip_files: A list of .zip files to decompress.
        output_dir: Optionally specify an output directory. Defaults to the present working directory.
        remove: Boolean value, which when set to True deletes level 1C .zip files after decompression is complete. Defaults to False.
    '''

    if type(zip_files) == str: zip_files = [zip_files]

    for zip_file in zip_files:
        assert zip_file[-4:] == '.zip', "Files to decompress must be .zip format."

    # Decompress each zip file
    for zip_file in zip_files:

        # Skip those files that have already been extracted
        if os.path.exists('%s/%s' % (output_dir, zip_file.split('\\')[-1].replace('.zip', '.SAFE'))):
            print('Skipping extraction of {}, as it has already been extracted in directory {}. If you want to re-extract it, delete the .SAFE file.'.format(zip_file, output_dir))

        else:
            print('Extracting {}'.format(zip_file))
            with zipfile.ZipFile(zip_file) as obj:
                obj.extractall(output_dir)
def RGBComposite(TileName,Tm,TileDir,OutDir):
    st = time.time()
    #L2A_201912_T49RGQ_20200426T222639_R10m_NIR
    R = glob.glob(os.path.join(TileDir,'L2A_{}_{}_*_R60m_B04.tif').format(Tm,TileName))[0]
    G = glob.glob(os.path.join(TileDir, 'L2A_{}_{}_*_R60m_B03.tif').format(Tm,TileName))[0]
    B = glob.glob(os.path.join(TileDir, 'L2A_{}_{}_*_R60m_B02.tif').format(Tm,TileName))[0]

    assert os.path.exists(R),"no file name {} exists".format(R)
    assert os.path.exists(G), "no file name {} exists".format(G)
    assert os.path.exists(B), "no file name {} exists".format(B)


    Rr = gdal.Open(R)
    Gr = gdal.Open(G)
    Br = gdal.Open(B)
    proj = Rr.GetProjection()
    geotrans = Rr.GetGeoTransform()
    basename = os.path.basename(R).split("_")
    Name = basename[0]+'_'+basename[1]+'_'+basename[2]+'_'+basename[4]+'_'+'RGB'+'.tif'
    Path = os.path.join(OutDir,Name)
    Size = Rr.RasterXSize
    Data = np.zeros(shape=(Size,Size,3))
    Data[:,:,0] = Rr.ReadAsArray()*0.0255
    del Rr
    Data[:, :, 1] = Gr.ReadAsArray()*0.0255
    del Gr
    Data[:, :, 2] = Br.ReadAsArray()*0.0255
    del Br

    Data = Image.fromarray(np.uint8(Data))
    Data =ImageOps.equalize(Data)

    write_Img(Data, Path, proj, geotrans, 1830, 1830, im_bands=3, dtype=gdal.GDT_Byte)
    del Data
    end = time.time()
    print(end - st)






InDir = r"D:\数据上线准备"
OutDir =r"D:\L2AV"
tileoutdir = r"D:\0RGB"
L2Alist = glob.glob(os.path.join(InDir,"*.zip"))

for L2A in L2Alist:
    decompress(L2A,OutDir)
    sets = glob.glob(os.path.join(OutDir,"*","*","*.tif"))
    TileName = sets[0].split("\\")[-2]
    Tm = os.path.basename(sets[0]).split("_")[1]
    Res = 60
    TileDir = glob.glob(os.path.join(OutDir,"*","*"))[0]
    RGBComposite(TileName,Tm,TileDir,tileoutdir)



    shutil.rmtree(glob.glob(os.path.join(OutDir,"*")))


