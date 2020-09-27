import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import gdal,osr,ogr
from scipy import signal
path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50_Anom.csv"


def ExtractVIFromCSV(path,SeasonType,var1,var2,year):

    assert SeasonType in ["Short","Long"], "SeasonType must be Short/Long"
    data = pd.read_csv(path)
    key1 = var1 + str(year) + SeasonType
    key2 = var2 + str(year) + SeasonType
    key3 = "PVI" + str(year) + SeasonType
    Var1List = data[key1].to_numpy()
    Var2List = data[key2].to_numpy()
    PVIList =  data[key3].to_numpy()
    COL = data["COL"].to_numpy()
    ROW = data["ROW"].to_numpy()

    return Var1List,Var2List,PVIList,COL,ROW
def GenValidRC():
    ref_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropMask\agg_clip50.tif"
    ref_raster = gdal.Open(ref_path)
    geo_t = ref_raster.GetGeoTransform()

    daShapefile = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropSHP\CropLandR.shp"
    # daShapefile = r"D:\Cornell\EthiopianDrought\westhihland\westhighland.shp"
    # daShapefile =r"D:\Cornell\EthiopianDrought\northhighlandshp\NorthHighland.shp"
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(daShapefile, 0)
    layer = dataSource.GetLayer()
    feature = layer.GetFeature(0)
    geo = feature.GetGeometryRef()
    geo = str(geo).split("((")[1].split("))")[0].split(",")
    x2 = []
    y2 = []
    for term in geo:
        x2.append(float(term.split(" ")[0]))
        y2.append(float(term.split(" ")[1]))

    x2 = np.array(x2)
    y2 = np.array(y2)
    x2 = (x2 - geo_t[0]) / geo_t[1]
    y2 = (y2 - geo_t[3]) / geo_t[5]
    x2 = x2.tolist()
    y2 = y2.tolist()

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id, p in enumerate(x2):
        ring.AddPoint(x2[id], y2[id])
    ring.AddPoint(x2[0], y2[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50_Anom.csv"
    Var1List, Var2List, Var3List,  COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", 2009)

    ROW2 = ROW.tolist()
    COL2 = COL.tolist()
    valid = []

    for id, r in enumerate(ROW2):

        point2 = ogr.Geometry(ogr.wkbPoint)
        point2.AddPoint(COL2[id], ROW2[id])

        if poly.Contains(point2):

            valid.append(id)
    return valid


valid = GenValidRC()
fig = plt.figure(figsize=(10, 6))
plt.xticks([])
plt.yticks([])
plt.axis('off')
SP = []
SPVI = []
LP = []
LPVI = []
SSIF=[]
LSIF=[]


for year in range(2007,2019):

    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", year)
    SP.append(Var2List[valid].mean())
    SPVI.append(PVIList[valid].mean())
    SSIF.append(Var1List[valid].mean())
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", year)
    LP.append(Var2List[valid].mean())
    LPVI.append(PVIList[valid].mean())
    LSIF.append(Var1List[valid].mean())
#
#
# SP = np.array(SP)
# SPVI = np.array(SPVI)
# SSIF = np.array(SSIF)
# LP = np.array(LP)
# LPVI = np.array(LPVI)
# LSIF = np.array(LSIF)
#
# #
# SP = (SP - SP.mean())/SP.std()
# SPVI = (SPVI - SPVI.mean())/SPVI.std()
# # SSIF = np.array(list(signal.detrend(SSIF)))
# SSIF = (SSIF - SSIF.mean())/SSIF.std()
# LP = (LP - LP.mean())/LP.std()
# LPVI = (LPVI - LPVI.mean())/LPVI.std()
# # LSIF = np.array(list(signal.detrend(LSIF)))
# LSIF = (LSIF - LSIF.mean())/LSIF.std()


ax = fig.add_subplot(2, 2, 1)
ax.set_xlabel("SIF Anomaly in Bega")
ax.set_ylabel("P Anomaly in Bega")
slope, intercept, r_value, p_value, std_err = stats.linregress(SSIF, SP)
print("SIF P",slope)
ax.scatter(SSIF,SP)
ax.text(0.2, 0.9,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
for i in range(12):
    ax.text(SSIF[i],SP[i],str(i+2007)[2:])
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xlabel("SIF Anomaly in bega")
ax2.set_ylabel("PVI Anomaly in bega")
slope, intercept, r_value, p_value, std_err = stats.linregress(SSIF, SPVI)
print("SIF PVI",slope)
ax2.scatter(SSIF,SPVI)
ax2.text(0.8, 0.9,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)
for i in range(12):
    ax2.text(SSIF[i],SPVI[i],str(i+2007)[2:])
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlabel("SIF Anomaly in krmet")
ax3.set_ylabel("P Anomaly in krmet")
slope, intercept, r_value, p_value, std_err = stats.linregress(LSIF, LP)
print("SIF P",slope)
ax3.scatter(LSIF,LP)
ax3.text(0.8, 0.2,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax3.transAxes)
for i in range(12):
    ax3.text(LSIF[i],LP[i],str(i+2007)[2:])
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xlabel("SIF Anomaly in krmet")
ax4.set_ylabel("PVI Anomaly in krmet")
slope, intercept, r_value, p_value, std_err = stats.linregress(LSIF, LPVI)
print("SIF PVI",slope)
ax4.scatter(LSIF,LPVI)
ax4.text(0.8, 0.2,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax4.transAxes)
for i in range(12):
    ax4.text(LSIF[i],LPVI[i],str(i+2007)[2:])
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_3.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()
