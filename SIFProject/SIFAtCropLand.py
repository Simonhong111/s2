import os
import glob
import numpy as np


def GetSIFAtCropLand(FilterCroplandGridPath,AllSIFPath,SaveSIFAtCropLandPath):

    Cropland = np.load(FilterCroplandGridPath)
    print("dd",Cropland["gridId"])
    SIF = np.load(AllSIFPath)
    print(Cropland.files)
    print(SIF.files)
    Date = []
    SIF757 =[]
    SIF771 = []
    SIFGridID = []
    Vec_SIF_AVG = []
    Vec_SIF_AVG_Norm = []
    EtRing = []
    GridTile = []
    waterP = []

    for idx,gridId in enumerate(SIF['sifgridId']):
        # print("gridId",gridId)
        for index2,gId in enumerate(Cropland["gridId"]):
            if int(gridId) == int(gId):
                print("gId",gId)
                Date.append(SIF["date"][idx])
                SIF757.append(SIF["sif757"][idx])
                SIF771.append(SIF["sif771"][idx])
                SIFGridID.append(SIF["sifgridId"][idx])
                Vec_SIF_AVG.append(SIF["vec_SIF_avg"][idx])
                Vec_SIF_AVG_Norm.append(SIF["vec_SIF_avg_norm"][idx])
                EtRing.append(SIF["EtRing"][idx])
                GridTile.append(Cropland["gridTile"][index2])
                waterP.append(Cropland["waterP"][index2])


    kwargs = {"date": Date, "sif757": SIF757, "sif771": SIF771, "sifgridId": SIFGridID,
              "vec_SIF_avg": Vec_SIF_AVG, "vec_SIF_avg_norm": Vec_SIF_AVG_Norm, \
              "EtRing": EtRing,'gridTile':GridTile, 'waterP':waterP}

    # np.savez(os.path.join(save_sif_path,str(dt.year)+"-"+str(dt.month))+"SifAgg.npz", **kwargs)
    np.savez(SaveSIFAtCropLandPath, **kwargs)



FilterCroplandGridPath = r"D:\Satellive\NPZ\MaskedGridMore.npz"
AllSIFPath = r"D:\Satellive\NPZ\2015-2018-SifAgg.npz"
SaveSIFAtCropLandPath =r"D:\Satellive\NPZ\2015-2018-SifAggAtCropLandMore.npz"  # 0.7 means the ratio of cropland is > 0.7
# GetSIFAtCropLand(FilterCroplandGridPath,AllSIFPath,SaveSIFAtCropLandPath)
ss = np.load(SaveSIFAtCropLandPath)
print(ss["date"])
print(len(ss["date"]))
i = 0
for date in ss["date"]:
    if "2018" in date:
        i += 1
print(i)
print(ss.files)

print(ss["gridTile"])