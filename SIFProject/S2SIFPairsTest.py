import numpy as np

from S2SIFPairs import *


s2point_path = r"D:\Sen2Projecton\Sen2Agg.npz"
sifpoint_path = r"D:\Sen2Projecton\SIF\sifnpz\2018-5.npz"
hb_shp_path = r"D:\Sen2Projecton\HBSHP\Hubei.shp"


Pairs(s2point_path,sifpoint_path,hb_shp_path,r"D:\Sen2Projecton\pairs.npz")

file = np.load(sifpoint_path)
print(file.files)

sif757 = file["sif757"]
for idx,itm in enumerate(sif757):
    print(file["date"][idx],file["sifgridId"][idx],file["sif757"][idx])
# gridId = file['gridId']
# gridEtRing = file['gridEtRing']
#
# date = file['Date']
# mean = file['mean']

# for idx , item in enumerate(gridId):
#     print(item,gridEtRing[idx],date[idx],mean[idx])

