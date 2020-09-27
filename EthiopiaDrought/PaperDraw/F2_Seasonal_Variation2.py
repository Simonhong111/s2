import numpy as np
import pandas as pd
from osgeo import gdal,osr,ogr
from matplotlib import pyplot as plt
from scipy import stats
path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"

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

    path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", 2009)

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

# ref_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropMask\agg_clip50.tif"
# ref_raster = gdal.Open(ref_path)
# plt.imshow(ref_raster.ReadAsArray())
# plt.scatter(COL[valid],ROW[valid])
# plt.show()
# for i in valid:
#     print(i,COL[i],ROW[i])
# breakpoint()


Year = 2009

# 1
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
plt.scatter(Var1List[valid],Var2List[valid])
plt.text(0.05, 0.95,'(a)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.text(1, 0.1,'R = {} P = {}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2a.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)

# 2
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
plt.scatter(Var1List[valid],PVIList[valid])
plt.text(1, 0.1,'R = {} P = {}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2b.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)

# 3
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
plt.scatter(Var1List[valid],Var2List[valid])
plt.text(1, 0.1,'R = {} P = {}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2c.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)

# 4
fig = plt.figure(figsize=(5, 3))
plt.xticks([])
plt.yticks([])
plt.text(0.05, 0.95,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
plt.scatter(Var1List[valid],PVIList[valid])
plt.text(1, 0.9,'R = {} P = {}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2d.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)

SP = []
SPV = []

SPVI = []
SPVIV = []

LP = []
LPV =[]

LPVI = []
LPVIV = []

for year in range(2007,2019):

    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
    SP.append(r_value)
    SPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
    SPVI.append(r_value)
    SPVIV.append(p_value)

    # long
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
    LP.append(r_value)
    LPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
    LPVI.append(r_value)
    LPVIV.append(p_value)

fig = plt.figure(figsize=(5, 3))
# plt.xticks([])
# plt.yticks([])

x = np.arange(12)
plt.bar(x-0.2,SP,width=0.4,label="P")
plt.bar(x+0.2,SPVI,width=0.4,label="PVI")
plt.legend()
plt.text(0.05, 0.95,'(e)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.xlabel('Year')
plt.ylabel("R")
plt.xticks(range(0,12,3),[str(2007+y)[2:] for y in range(0,12,3)])
plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2e.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)

fig = plt.figure(figsize=(5, 3))
# plt.xticks([])
# plt.yticks([])

x = np.arange(12)
plt.bar(x-0.2,LP,width=0.4,label="P")
plt.bar(x+0.2,LPVI,width=0.4,label="PVI")
plt.legend()
plt.text(0.05, 0.95,'(f)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = plt.gca().transAxes)
plt.xlabel('Year')
plt.ylabel("R")
plt.xticks(range(0,12,3),[str(2007+y)[2:] for y in range(0,12,3)])

plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2f.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)


# plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
# plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()
