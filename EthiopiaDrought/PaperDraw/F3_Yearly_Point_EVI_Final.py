import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import gdal,osr,ogr
from scipy import signal
import matplotlib.gridspec as gridspec
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

SP = []
SPVI = []
LP = []
LPVI = []
SSIF=[]
LSIF=[]


for year in range(2007,2019):

    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "EVI", "RF", year)
    SP.append(Var2List[valid].mean())
    SPVI.append(PVIList[valid].mean())
    SSIF.append(Var1List[valid].mean())
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "EVI", "RF", year)
    LP.append(Var2List[valid].mean())
    LPVI.append(PVIList[valid].mean())
    LSIF.append(Var1List[valid].mean())

myear = [y for y in range(2007,2019)]
print("*sp",list(zip(myear,SP)))
print("*spvi",list(zip(myear,SPVI)))
print("*ssif",list(zip(myear,SSIF)))
print("*lp",list(zip(myear,LP)))
print("*lpvi",list(zip(myear,LPVI)))
print("*lsif",list(zip(myear,LSIF)))

fig2 = plt.figure(constrained_layout=True,figsize=(10,10))
spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig2)
f2_ax1 = fig2.add_subplot(spec2[0, 0])
f2_ax2 = fig2.add_subplot(spec2[0, 1])
f2_ax3 = fig2.add_subplot(spec2[1, 0])
f2_ax4 = fig2.add_subplot(spec2[1, 1])


f2_ax1.set_xlabel("EVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax1.set_ylabel("P Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
slope, intercept, r_value, p_value, std_err = stats.linregress(SSIF, SP)
print("SIF P",slope)
f2_ax1.scatter(SSIF,SP)
f2_ax1.text(0.2, 0.9,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax1.transAxes)
f2_ax1.text(0.5, 0.95,'(a)'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax1.transAxes)
for i in range(12):
    f2_ax1.text(SSIF[i],SP[i],str(i+2007)[2:])


f2_ax2.set_xlabel("EVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax2.set_ylabel("PVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
slope, intercept, r_value, p_value, std_err = stats.linregress(SSIF, SPVI)
print("SIF PVI",slope)
f2_ax2.scatter(SSIF,SPVI)
f2_ax2.text(0.8, 0.9,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax2.transAxes)
f2_ax2.text(0.5, 0.95,'(b)'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax2.transAxes)
for i in range(12):
    f2_ax2.text(SSIF[i],SPVI[i],str(i+2007)[2:])

f2_ax3.set_xlabel("EVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax3.set_ylabel("P Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
slope, intercept, r_value, p_value, std_err = stats.linregress(LSIF, LP)
print("SIF P",slope)
f2_ax3.scatter(LSIF,LP)
f2_ax3.text(0.8, 0.2,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax3.transAxes)
f2_ax3.text(0.5, 0.95,'(c)'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax3.transAxes)
for i in range(12):
    f2_ax3.text(LSIF[i],LP[i],str(i+2007)[2:])

f2_ax4.set_xlabel("EVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax4.set_ylabel("PVI Anomaly",fontdict={'fontname':'Times New Roman','fontsize':16})
slope, intercept, r_value, p_value, std_err = stats.linregress(LSIF, LPVI)
print("SIF PVI",slope)
f2_ax4.scatter(LSIF,LPVI)
f2_ax4.text(0.8, 0.2,'R = {} P = {}'.format(round(r_value, 3),round(p_value, 3)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax4.transAxes)
f2_ax4.text(0.5, 0.95,'(d)'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax4.transAxes)
for i in range(12):
    f2_ax4.text(LSIF[i],LPVI[i],str(i+2007)[2:])
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_3EVI.png',bbox_inches='tight',dpi=fig2.dpi,pad_inches=0.05)
plt.show()

