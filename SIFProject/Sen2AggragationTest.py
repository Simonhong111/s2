from Sen2Aggragation import *
from Sen2Clip import *
import numpy as np
grid_path = r"D:\Sen2Projecton\NPZ\MaskedGrid.npz"

sen2dir = r"J:\2018S2docker"

save_npz_path  = r"D:\Sen2Projecton\Sen2Agg.npz"

SaveSen2Agg(grid_path,sen2dir,save_npz_path)




# sen2agg = np.load(save_npz_path)
#
# print(sen2agg.files)
#
#
#
# gridId = sen2agg['gridId']
# gridEtRing = sen2agg[ 'gridEtRing']
# print(len(gridEtRing))
# Date = sen2agg[ 'Date']
# mean = sen2agg[ 'mean']
# EtRing = sen2agg["gridEtRing"]
# s2path = sen2agg['s2path']
# print(s2path.__len__())
# print(gridId.__len__())
#
# for idx,item in enumerate(gridId):
#     mSenL2A = SenL2AReader(s2path[idx])
#     band = mSenL2A.getImgData("B03")
#
# #
#     ulx,uly = EtRing[idx][0],EtRing[idx][1]
#     lrx,lry = EtRing[idx][4],EtRing[idx][5]
#
#
#
#
#     del band
#     respath = os.path.join(mSenL2A.getImgDataPath(), "R" + str(20) + "m")
#
#     # 定位到具体波段，如 B02 波段
#     rasterpath = glob.glob(os.path.join(respath, '*' + "B12" + "*"))
#     qipath = glob.glob(os.path.join(mSenL2A.getQiDataPath(), '*' + "MSK_CLDPRB" + "*"))
#     GetRoi(rasterpath[0],r"D:\Sen2Projecton\Test\test.tif",[uly,ulx,lry, lrx])
#     GetRoi(qipath[0],r"D:\Sen2Projecton\Test\testcld.tif",[uly,ulx,lry, lrx])
#     raster = gdal.Open(r"D:\Sen2Projecton\Test\test.tif",0).ReadAsArray()
#     cld =  gdal.Open(r"D:\Sen2Projecton\Test\testcld.tif",0).ReadAsArray()
#     mask = (raster < 10000) & (raster > 0) & (cld < 30)
#     mask = np.where(mask == True)
#     mMean = np.mean(raster[mask])*0.0001
#     print("mean",mMean,mean[idx],s2path[idx])








